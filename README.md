# dicom_to_nrrd
Simple program for converting dicom images and structure sets to nrrd files. This will create mask images of structures corresponding to an image series, and exported nrrd files can then be used for pyradiomics. 

Clone this repository, or download the 4 python files and place them in one directory. Run dicom_to_nrrd.py for the program.

install all dependencies in requirements.txt (ex. pip install -r requirements.txt)

The program will prompt you for directories to find files and where to save files, but it is best to provide these as arguments when running the program. 
There are 5 arguments (2 needed minimum):
--img_dir: specify the folder containing a single patient's image series.
--structure_dir: specify the folder containing a single patient's RTSTRUCT file.
--modality: specify image modality you wish to convert (Default is CT. Other options are PET or MR)
--save_dir: specify where you wish to save your nrrd files. if You wish to place converted structures in a separate folder, also specify the argument "save_dir_struct\"
--save_dir_struct: Specify a separate folder you'd like to save the struct nrrd files to.

Example:

on the command line (all one line) 
>> python dicom_to_nrrd.py --img_dir path/to/dicom/img/series --structure_dir path/to/dicom/structure/files  --modality MR --save_dir path/to/save/nrrd/img --save_dir_struct path/to/save/nrrd/structure

Note that it is ok to have the structure files in the same directory as the images. In this case you only need to include the image directory, img_dir. Similarly, it is ok to save all data to a single directory, and you then only need to include --save_dir.

Example 2:

to incorporate this into your python program, you just need to call the dicom_to_nrrd.convert_all_dicoms_to_nrrd(img_dir, struct_dir, modality, save_paths=None) function.

For now, you must have only one patient's data in the directories, and also only one image series. 
