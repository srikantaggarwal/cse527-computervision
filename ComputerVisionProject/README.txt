The code is implemented by Arvind Chaudahry, Srikant Aggarwal, Roman Pavlushchenko for CSE527 course project.
Requirements:
	1. Install quadprod optimization tool.
	2. Add  directory ImageRecH_V01 to the path.
For RGB images call RemoveShadowRGB and for gray-scale images call RemoveShadowGray.
For grayscale images, provide grayscale images only.
Main function:
	RemoveShadow(InputImage, MaskImage)
		This function reads InputImage that contains the shadow and MaskImage that is corresponding mask for the image.
		Now it calculates the optimal and final illumination change model for the penumbra area. Using quadprog.
		Now it calls RemoveShadowEffect.
		After removing the shadow effect, now we call extendedMask and calcGradientST to reconstruct the texture in the gradient field in penumbra area. Once we have the gradients, we reconstruct the shadow free image by using Poisson solver.
	RemoveShadowEffect(modelsX,  modelsY, Gx, Gy, R)
		This function removes the shadow effect from Gradients Gx and Gy using the illumination change models calculated in previous function.
		
External code:
	quadprog
	poisson solver.
	bwtraceboundary

NOTE: For some images, bwtraceboundary function doesn't work.
We have attached the results as well in two directories named 'Colored' and 'grayscale'.
	