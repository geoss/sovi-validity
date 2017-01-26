import numpy as np
from scipy.stats.mstats import zscore as ZSCORE
import mdp as MDP  # causing sklearn deprication warning
from operator import itemgetter
import copy


class SPSS_PCA:
	'''
	A class that integrates most (all?) of the assumptions SPSS imbeds in their
    implimnetation of principle components analysis (PCA), which can be found in
    thier GUI under Analyze > Dimension Reduction > Factor. This class is not
	intended to be a full blown recreation of the SPSS Factor Analysis GUI, but
	it does replicate (possibly) the most common use cases. Note that this class
	will not produce exactly the same results as SPSS, probably due to differences
	in how eigenvectors/eigenvalues and/or singular values are computed. However,
	this class does seem to get all the signs to match, which is not really necessary
	but kinda nice. Most of the approach came from the official SPSS documentation.

	References
	----------
	ftp://public.dhe.ibm.com/software/analytics/spss/documentation/statistics/20.0/en/client/Manuals/IBM_SPSS_Statistics_Algorithms.pdf
	http://spssx-discussion.1045642.n5.nabble.com/Interpretation-of-PCA-td1074350.html
	http://mdp-toolkit.sourceforge.net/api/mdp.nodes.WhiteningNode-class.html
	https://github.com/mdp-toolkit/mdp-toolkit/blob/master/mdp/nodes/pca_nodes.py

	Parameters
	----------
	inputs:  numpy array
			 n x k numpy array; n observations and k variables on each observation
	reduce:  boolean (default=False)
			 If True, then use eigenvalues to determine which factors to keep; all
			 results will be based on just these factors. If False use all factors.
	min_eig: float (default=1.0)
			 If reduce=True, then keep all factors with an eigenvalue greater than
			 min_eig. SPSS default is 1.0. If reduce=False, then min_eig is ignored.
	varimax: boolean (default=False)
			 If True, then apply a varimax rotation to the results. If False, then
			 return the unrotated results only.

	Attributes
	----------
	z_inputs:	numpy array
				z-scores of the input array.
	comp_mat:	numpy array
				Component matrix (a.k.a, "loadings").
	scores:		numpy array
				New uncorrelated vectors associated with each observation.
	eigenvals_all:	numpy array
				Eigenvalues associated with each factor.
	eigenvals:	numpy array
				Subset of eigenvalues_all reflecting only those that meet the
				criterion defined by parameters reduce and min_eig.
	weights:    numpy array
				Values applied to the input data (after z-scores) to get the PCA
				scores. "Component score coefficient matrix" in SPSS or
				"projection matrix" in the MDP library.
	comms: 		numpy array
				Communalities
	sum_sq_load: numpy array
				 Sum of squared loadings.
	comp_mat_rot: numpy array or None
				  Component matrix after rotation. Ordered from highest to lowest
				  variance explained based on sum_sq_load_rot. None if varimax=False.
	scores_rot:	numpy array or None
				Uncorrelated vectors associated with each observation, after
				rotation. None if varimax=False.
	weights_rot: numpy array or None
				Rotated values applied to the input data (after z-scores) to get
				the PCA	scores. None if varimax=False.
	sum_sq_load_rot: numpy array or None
				 Sum of squared loadings for rotated results. None if
				 varimax=False.

	'''

	def __init__(self, inputs, reduce=False, min_eig=1.0, varimax=False):
		z_inputs = ZSCORE(inputs)  # seems necessary for SPSS "correlation matrix" setting (their default)

		# run base SPSS-style PCA to get all eigenvalues
		pca_node = MDP.nodes.WhiteningNode()  # settings for the PCA
		scores = pca_node.execute(z_inputs)  # base run PCA
		eigenvalues_all = pca_node.d   # rename PCA results

		# run SPSS-style PCA based on user settings
		pca_node = MDP.nodes.WhiteningNode(reduce=reduce, var_abs=min_eig)  # settings for the PCA
		scores = pca_node.execute(z_inputs)  # run PCA  (these have mean=0, std_dev=1)
		weights = pca_node.v  # rename PCA results (these might be a transformation of the eigenvectors)
		eigenvalues = pca_node.d   # rename PCA results
		component_matrix = weights * eigenvalues  # compute the loadings
		component_matrix = self._reflect(component_matrix)   # get signs to match SPSS
		communalities = (component_matrix**2).sum(1)   # compute the communalities
		sum_sq_loadings = (component_matrix**2).sum(0) # note that this is the same as eigenvalues
		weights_reflected = component_matrix/eigenvalues  # get signs to match SPSS
		scores_reflected = np.dot(z_inputs, weights_reflected)  # note that abs(scores)=abs(scores_reflected)

		if varimax:
			# SPSS-style varimax rotation prep
			c_normalizer = 1. / MDP.numx.sqrt(communalities)  # used to normalize inputs to varimax
			c_normalizer.shape = (component_matrix.shape[0],1)  # reshape to vectorize normalization
			cm_normalized = c_normalizer * component_matrix  # normalize component matrix for varimax

			# varimax rotation
			cm_normalized_varimax = self._varimax(cm_normalized)  # run varimax
			c_normalizer2 = MDP.numx.sqrt(communalities)  # used to denormalize varimax output
			c_normalizer2.shape = (component_matrix.shape[0],1)  # reshape to vectorize denormalization
			cm_varimax = c_normalizer2 * cm_normalized_varimax  # denormalize varimax output

			# reorder varimax component matrix
			sorter = (cm_varimax**2).sum(0)  # base the ordering on sum of squared loadings
			sorter = zip(sorter.tolist(), range(sorter.shape[0]))  # add index to denote current order
			sorter = sorted(sorter, key=itemgetter(0), reverse=True)  # sort from largest to smallest
			sum_sq_loadings_varimax, reorderer = zip(*sorter)  # unzip the sorted list
			sum_sq_loadings_varimax = np.array(sum_sq_loadings_varimax)  # convert to array
			cm_varimax = cm_varimax[:,reorderer]  # reorder component matrix

			# varimax scores
			cm_varimax_reflected = self._reflect(cm_varimax)  # get signs to match SPSS
			varimax_weights = np.dot(cm_varimax_reflected,
							  np.linalg.inv(np.dot(cm_varimax_reflected.T,
							  cm_varimax_reflected))) # CM(CM'CM)^-1
			scores_varimax = np.dot(z_inputs, varimax_weights)
		else:
			comp_mat_rot = None
			scores_rot = None
			weights_rot = None

		# assign output variables
		self.z_inputs = z_inputs
		self.scores = scores_reflected
		self.comp_mat = component_matrix
		self.eigenvals_all = eigenvalues_all
		self.eigenvals = eigenvalues
		self.weights = weights_reflected
		self.comms = communalities
		self.sum_sq_load = sum_sq_loadings
		self.comp_mat_rot = cm_varimax_reflected
		self.scores_rot = scores_varimax
		self.weights_rot = varimax_weights
		self.sum_sq_load_rot = sum_sq_loadings_varimax

	def _reflect(self, cm):
		# reflect factors with negative sums; SPSS default
		cm = copy.deepcopy(cm)
		reflector = cm.sum(0)
		for column, measure in enumerate(reflector):
			if measure < 0:
				cm[:,column] = -cm[:,column]
		return cm

	def _varimax(self, Phi, gamma = 1.0, q = 100, tol = 1e-6):
		# downloaded from http://en.wikipedia.org/wiki/Talk%3aVarimax_rotation
		# also here http://stackoverflow.com/questions/17628589/perform-varimax-rotation-in-python-using-numpy
		p,k = Phi.shape
		R = np.eye(k)
		d=0
		for i in range(q):
			d_old = d
			Lambda = np.dot(Phi, R)
			u,s,vh = np.linalg.svd(np.dot(Phi.T,np.asarray(Lambda)**3 - (gamma/p) *
							np.dot(Lambda, np.diag(np.diag(np.dot(Lambda.T,Lambda))))))
			R = np.dot(u,vh)
			d = np.sum(s)
			if d_old!=0 and d/d_old < 1 + tol:
				break
		return np.dot(Phi, R)
