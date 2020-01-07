
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from skimage import data, img_as_float
import matplotlib.image as mpimg
from skimage.segmentation import (morphological_chan_vese,
                                  morphological_geodesic_active_contour,
                                  inverse_gaussian_gradient,
                                  checkerboard_level_set)
from PIL import Image, ImageEnhance
#import mat4py
import scipy.io as io
import xlrd
import skimage
import PIL

def store_evolution_in(lst):
    """Returns a callback function to store the evolution of the level sets in
    the given list.
    """

    def _store(x):
        lst.append(np.copy(x))

    return _store


# Morphological ACWE
#image = img_as_float(data.camera())
#image1=mpimg.imread('Tumor1.png')
def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.299, 0.587, 0.114])

def active_contour(image, snake, alpha=0.01, beta=0.1,
                   w_line=0, w_edge=1, gamma=0.01,
                   bc='periodic', max_px_move=1.0,
                   max_iterations=2500, convergence=0.1):
    return

'''

#img = data.astronaut()
#img = rgb2gray(img)
img1 = mpimg.imread('tum4.png')
image = rgb2gray(img1)
plt.figure()
plt.imshow(image, cmap = plt.get_cmap('gray'))
image= change_contrast(Image.open('tum4.png'), 50)

#plt.figure()
plt.imshow(x)
plt.show()
'''


for a in range (5,6):

    im = Image.open('tum'+str(a)+'.png')
    enhancer = ImageEnhance.Brightness(im)
    enhanced_im = enhancer.enhance(0.5)
    enhancer = ImageEnhance.Contrast(enhanced_im)
    enhanced = enhancer.enhance(0.2)
    image1=enhanced.convert('L')
    image = np.asarray(image1)
    #plt.figure()
    #plt.imshow(enhanced)

    #enhanced_im.save("enhanced.sample5.png")

    #plt.show()
    #image = img_as_float(image1)

    # Initial level set
    init_ls = checkerboard_level_set(image.shape, 15)
    # List with intermediate results for plotting the evolution
    evolution = []
    callback = store_evolution_in(evolution)
    ls = morphological_chan_vese(image, 140, init_level_set=init_ls, smoothing=10, lambda1=1.2, lambda2=1.5,
                                 iter_callback=callback)

    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    ax = axes.flatten()


    ax[0].imshow(image,    cmap="gray")
    ax[0].set_axis_off()
    ax[0].contour(ls, [0.8], colors='b')
    ax[0].set_title("Morphological ACWE segmentation", fontsize=12)



    
    ax[1].imshow(ls, cmap="gray")
    ax[1].set_axis_off()
    contour = ax[1].contour(evolution[15], [0.5], colors='g')
    contour.collections[0].set_label("Iteration 2")
    contour = ax[1].contour(evolution[7], [0.5], colors='y')
    contour.collections[0].set_label("Iteration 7")
    contour = ax[1].contour(evolution[-1], [0.5], colors='r')
    contour.collections[0].set_label("Iteration 35")
    ax[1].legend(loc="upper right")
    title = "Morphological ACWE evolution"
    ax[1].set_title(title, fontsize=12)
    ls1=np.asarray(ls)
    ls1 = 255*ls1
    img3 = PIL.Image.fromarray(ls1)
    
    print(len(ls1))
    print(type(img3))
    plt.figure()
    plt.imshow(img3)

    img = Image.new('RGB', (512, 512),'red')
    image_copy = img.copy()
    image_copy.paste(img3, (20, 20))
    plt.figure()
    plt.imshow(image_copy)


    # Morphological GAC
    #image = img_as_float(data.coins())
    gimage = inverse_gaussian_gradient(image)

    # Initial level set
    init_ls = np.zeros(image.shape, dtype=np.int8)
    init_ls[10:-10, 10:-10] = 1

    # List with intermediate results for plotting the evolution
    evolution = []
    callback = store_evolution_in(evolution)
    ls = morphological_geodesic_active_contour(gimage, 230, init_ls,
                                               smoothing=1, balloon=-1,
                                               threshold=0.5,
                                               iter_callback=callback)

    ax[2].imshow(image, cmap="gray")
    ax[2].set_axis_off()
    ax[2].contour(ls, [0.5], colors='r')
    ax[2].set_title("Morphological GAC segmentation", fontsize=12)

    io.savemat('temp',{"foo":ls})

    ax[3].imshow(ls, cmap="gray")
    ax[3].set_axis_off()
    contour = ax[3].contour(evolution[0], [0.5], colors='g')
    contour.collections[0].set_label("Iteration 0")
    contour = ax[3].contour(evolution[100], [0.5], colors='y')
    contour.collections[0].set_label("Iteration 100")
    contour = ax[3].contour(evolution[-1], [0.5], colors='r')
    contour.collections[0].set_label("Iteration 230")
    ax[3].legend(loc="upper right")
    title = "Morphological GAC evolution"
    ax[3].set_title(title, fontsize=12)
    fig.tight_layout()




plt.show()
