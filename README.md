# IDL-codes
These are the IDL codes used in my position as a research assistant.

star_calib_millstone_tim1.pro and unwarp_millstone_tim1.pro are both codes for which I started with existing code to perform a certain tasks, and it was up to me to improve them by making them more easily readable, easier to edit, and adding new features to them.

star_calib_millstone_tim1.pro is used for calibration. Initially, images being used are highly distorted with a fish-eye like effect. This code searches for locations of stars in the image; We can easily calculate where in the sky they should be, so once we find where in the image they are, we can create a distortion function that applies to all images taken with this camera. It also determines how sensitive the imager is to brightness. Again, we know how bright the stars are supposed to be, so we can determine how much brightness each pixel value corresponds to.

unwarp_millstone_tim1.pro uses the calibration data from star_calib_millstone_tim1.pro. It uses the distortion function to unwarp the images so that each pixel represents a constant latitude/longitude interval. It also adjusts the pixel values such that the pixel value is equal to the value in the brightness unit Rayleighs. It also has capabilites to create a plot known as a "velogram" that aligns cross sections of each image unwarped in a given data set to determine the speed and direction of propagation of different events and features. It also has capabilities to plot the locations of satellites passing overhead, in order to compare data from our images to data taken by those satellites. 

Each routine calls numerous subroutines which are not included in this repository.
