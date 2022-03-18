# imageopt
A handy tool for VCNEB calculations: reads the Images file and (greedily) optimizes the atomic positions
Simply run optimize_images.py in the same folder where the Images file is.
The script reads the two POSCAR files and attempts to:
  1) Rotate the crystal axes of Image_end as close as possible to those of Image_ini
  2) Swap the coordinates of the atoms in Image_end to make them as close as possible to those of the same atomic specie in Image_ini
Then writes a Images_Order file.
