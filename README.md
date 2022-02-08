# Settlement Time Calc Module

## Purpose
 Using both Terzaghi's and Barron's analytical solutions, this script returns expected and residual settlement magnitudes given general soil properties and duration. In today's realm of computational methods which can be highly precise yet convoluted, analytical solutions offer good estimates that serves as a crosscheck. If calibrated properly, the results from this module can be comparable to numerically intensive methods (e.g. FDM, FEM)

## Input
 Inputs required from the user include strata thickness, compressibility parameters, stress states and many more. For more information, refer to the annotated Jupyter Notebook. If compressibility parameters are not available, the function is calibrated to use typical [Singapore Upper Marine Clay properties](https://www.sciencedirect.com/science/article/pii/S003808061500061X)(Myint et al,2015)

## Output
 Several functions are included in this module which will be elaborated on in further updates. Nevertheless, the main function is to calculate the expected and residual settlement within a given duration. 
By default, the script produces a Settlement vs Time curve. The Terzaghi and Barron charts used can also be displayed if needed.

## Future Development
1. Incorporate in multi-layer stratum and less common isochrone situations (e.g. non-uniform PWP with single way drainage)
2. Estimate clay strength improvement given original strength/SPT/CPT values.
