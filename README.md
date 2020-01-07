# Chan_Vese-algorithm-for-segmentattion
After extracting the tumor region from the original image, it is necessary to do a segmentation of the tumor and produce a tumor mask. As the first step of the segmentation process image preprocessing was done. In there, adjust the image contrast, brightness level and the sharpness of the images. After applying the preprocessing techniques Morphological Active Contours (Morphological chan verse) technique was applied to identify the tumor region. 

After segment the tumor, a tumor mask was generated. Then the generated mask is compared with original tumor mask by ising few paramertes.
comparison_masks.m file includes mask comparission and Contour Detection results.pdf file includes the results for several test cases. 
