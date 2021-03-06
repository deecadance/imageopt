#! /usr/bin/python

import numpy as np
from scipy.spatial.transform import Rotation
filesave = []

rotate = True
reorder = True
## LATTICE PARAMETER
lattpar_a = np.zeros((3,3))
lattpar_b = np.zeros((3,3))
############################################################
########  PART 1: IMPORT DATA FROM Images FILE  ############
############################################################
## READ THE Images FILE AND STORE IT FOR LATER
## ALSO: READ NUMBER OF TYPES OF ATOMS AND NUMBER OF ATOMS
with open('Images') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        filesave.append(line.strip().split())
        if i == 5:
            elements = line.strip().split()
            tot_elements = len(elements)
            print(elements)
        if i == 6:
            num_elements = np.asarray(line.strip().split()).astype("int")
            print(num_elements)

## ARRAY WITH COORDS OF Image_ini
coords_elem_a = np.zeros((np.sum(num_elements),3))
## ARRAY WITH THE TYPES OF ELEMENTS REPEATED, eg. H2O --> ["H", "H", "O"]
list_elem = []
for i, element in enumerate(elements):
   list_elem.append(np.repeat(element,num_elements[i]))
list_elem = np.concatenate(list_elem)
##print(list_elem)
## ARRAY WITH COORDS OF Image_end
coords_elem_b = np.zeros((np.sum(num_elements),3))

## READ THE Images AGAIN, NOW STORE Image_ini AND Image END INTO TWO SEPARATE ARRAYS
with open('Images') as f:
    lines = f.readlines()
    j, k = 0, 0
    ## FIRST SEVEN LINE IN POSCAR ARE RESERVED FOR HEADER AND LATTICE
    for i, line in enumerate(lines):
        if (i >= 2)  & (i <= 4):
            lattpar_a[:, i-2] = np.asarray((line.strip().split())).astype("float")
        if (i >= 8 ) & (i < 8+np.sum(num_elements)):
            coords_elem_a[j,:] =  np.asarray(line.strip().split()).astype("float")
            j = j + 1
        if (i >= 10+np.sum(num_elements))  & (i <= 12+np.sum(num_elements)):
            lattpar_b[:, i-(10+np.sum(num_elements))] = np.asarray((line.strip().split())).astype("float")
        if (i >= 16+np.sum(num_elements)) & (i < 16+2*np.sum(num_elements)):
            coords_elem_b[k,:] =  np.asarray(line.strip().split()).astype("float")
            k = k + 1
############################################################
########  PART 2: ROTATE LATTPAR_B -> LATTPAR_A  ###########
############################################################
## SWITCH TO CARTESIAN COORDINATES (NEEDED FOR ROTATION!)
if rotate:
    ## FIND BEST ROTATION TO BRING LATTPAR_B TO LATTPAR_A (see scipy.spatial.transform.Rotation.align_vectors)
    rotation, rmsd = Rotation.align_vectors(lattpar_a[:,:].T,lattpar_b[:,:].T)
    ## MATRIX PRODUCT: LATTPAR_B --> ROTATED TO BE CLOSE TO LATTPAR_A
    old_lattpar_b = lattpar_b
    lattpar_b = np.dot(rotation.as_matrix(),lattpar_b)

    cart_coords_elem_a = np.dot(coords_elem_a,lattpar_a)
    cart_coords_elem_b = np.dot(coords_elem_b,lattpar_b)
    coords_output_a = np.dot(cart_coords_elem_a,np.linalg.inv(lattpar_a))
    coords_output_b = coords_elem_b #coords_rotate_b 
    #print(coords_rotate_b)
############################################################
########  PART 3: OPTIMIZE COORDS OF POSCAR B  #############
############################################################
## LIST OF MASKS CORRESPONDING TO THE DIFFERENT ELEMENTS
if reorder:
    mask_elem_list = np.full((len(list_elem),tot_elements),True,dtype=bool)
    for i, element in enumerate(elements):
        mask_elem_list[:,i] = list_elem == element
    
    coords_output_b = np.ones((np.sum(num_elements),3))
    ## K KEEPS TRACK OF WHICH LINE WE'RE ON ON THE POSITIONS
    k = 0
    ## LOOP OVER ELEMENTS (USING THE MASK)
    for l in range(mask_elem_list.shape[1]):
        ## FILL MATRIX WITH DISTANCES BETWEEN ATOMS IN THE TWO POSCARS
        ## THE DISTANCE SHOULD USE CARTESIAN COORDINATES!
        dist = np.zeros((len(cart_coords_elem_a[mask_elem_list[:,l]]), len(cart_coords_elem_a[mask_elem_list[:,l]])))
        for i, line_a in enumerate(cart_coords_elem_a[mask_elem_list[:,l]]): #enumerate(coords_elem_a[mask_elem_list[:,l]]):
            for j, line_b in enumerate(cart_coords_elem_b[mask_elem_list[:,l]]): #enumerate(coords_elem_b[mask_elem_list[:,l]]):
                dist[i,j] = np.sqrt(np.sum((line_a[:]-line_b[:])**2))
#        print(dist)
        ## LOOP OVER DISTANCES. ALGORITHM IS GREEDY: PUTS THE CLOSEST ATOM AND NEVER CHECKS IT AGAIN
        ## ALL THE "new_unique_coord" MESS IS BECAUSE WE CANNOT USE THE SAME POSITION TWICE
        for i, dist_line in enumerate(dist):
    #        print(np.argsort(dist_line), dist_line)	
            ctrl = 0
    	## FLAG TO DETERMINE WHETHER THIS COORD WAS ALREADY USED
            new_unique_coord = False
            while (not new_unique_coord) and (ctrl < len(coords_elem_b[mask_elem_list[:,l]])):
                if np.any(np.all(coords_elem_b[mask_elem_list[:,l]][np.argsort(dist_line)][ctrl,:] == coords_output_b, axis=1)):
                    ctrl = ctrl + 1
                else:
                    coords_output_b[k,:] = coords_elem_b[mask_elem_list[:,l]][np.argsort(dist_line)][ctrl,:]
                    new_unique_coord = True
##HERE
                    print(dist_line[np.argsort(dist_line)][ctrl], dist_line, ctrl)
            ## IF AT THE END OF THE LOOP WE COULD NOT FIND A NEW COORDINATE, THERE IS SOMETHING WRONG
            if new_unique_coord == False:
                print("ERROR! CHECK YOUR COORDINATES!")
                sys.exit(1)
            k = k + 1

## WRITE THE NEW FILE ALREADY PARSED AND READY TO USE
with open('Images_Order', 'w') as f:
    j = 0
    for i, line in enumerate(filesave):
        if isinstance(line[0], (int, float)):
            f.write('       '.join(["{:.16f}".format(elem) for elem in line])+"\n")
        elif (i >= 2)  & (i <= 4):
            f.write('       '.join(["{:.16f}".format(elem) for elem in lattpar_a[:,i-2]])+"\n")
        elif (i >= 8 ) & (i < 8+np.sum(num_elements)):
            f.write('       '.join(["{:.16f}".format(elem) for elem in coords_output_a[i-8,:]])+"\n")
        elif (i >= 10+np.sum(num_elements))  & (i <= 12+np.sum(num_elements)):
            f.write('       '.join(["{:.16f}".format(elem) for elem in lattpar_b[:,i-(10+np.sum(num_elements))]])+"\n")
        elif (i >= 16+np.sum(num_elements)) & (i < 16+2*np.sum(num_elements)):
            f.write('       '.join(["{:.16f}".format(elem) for elem in coords_output_b[j,:]])+"\n")
            j = j + 1
        else:
            f.write('       '.join([str(elem) for elem in line])+"\n")
