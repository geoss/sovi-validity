# sovi-validity
Analysis of logical consistence of social vulnerability indices

This is a script to calculate the Social Vulnerability (Sovi) Index. The Sovi index measures a place's social vulnerability to natural hazards and has been used in hundreds of publications, in both the academic and policy domains.  Social vulnerability is the idea that different populations experience the same hazard event in different ways.  For example, Katrina showed that some socio-economic groups were more vulnerable than others.

The code here computes the SoVI Index at multiple geographic scales and using different sets of inputs. We do this in an effort to measure the [convergent](https://en.wikipedia.org/wiki/Convergent_validity) and [construct validity](https://en.wikipedia.org/wiki/Construct_validity) of the index.  

*Convergent Validity* is the idea that multiple measurements of the same thing, using a valid instruments, should yield similar absolute and relative measurements. For example two thermometers measuing the same cup of water should yield the same approximate temperature- this would be an example of validitiy of absolute measurement. An less rigorous concept of validity is to consider three cups ordered from hottest (A) to coldest (C), the two thermometers would be valid in a relative sense if their measurements of cups A, B, C differed absolutely (they measured different temperatures for each cup) but still placed cup A as the warmest and cup C as the coldest.

This code shows in this analysis that the Sovi Index fails this "cup test" that is it often mixes up the orders of the cup. Counties in the United States that are high vulnerability at one scale of measurement (or form the index) are often low vulnerability in a slightly different version of the index.

*Construct Validity* is the idea that a instrument correctly measures the thing it purports to measure.  Data from the thermometer in the previous example would lack construct validity if it measured air pressure instead of temperature.  The construct validity of an index is best interrogated by examining how variable combine to form the index.  That what a social index means can be understood to examining its ingredients.  We examine the construct validity of the SoVI index by looking at how variables contribute to the index.  Looking at the variable-by-variable contribution to the index produceses all sorts of nonsenecial insights: like wealthier places are more socially vulnerable.   

## Variables and Components
The Sovi Index is constructed using a tecnhnique called Principal Components Analysis, this is matrix decomposition method that uses the covariance matrix of the input data. Usually, in the social sciences one treats the "compents" what come of out a PCA as latent variables. For example, in Sovi it comon to fine components that measure things like "race and class". In this analysis we also show that this latent variable approach has maked some underlying problems with the Soci index, namely that variables contribute to the index in ways that are profoundly counter intuitive.

##There is a paper
For an in-depth discussion of these ideas please see the companion paper to this anlysis URL or contact the suthors.

##Data Prep
All of the code necessary to transform raw data from the ACS inthe 28 variables used to construct SoVI is included.  There is a lot of data wrangling in this code, combining variables into new variables, computing standard errors, etc. Its all farily straightforward and we will not walk through the details here.  The code should be self-explanatory.
