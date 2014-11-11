1) Run the executable. It will create a folder ./IO

2) Copy some images to the IO folder. Filenames MUST be of the form 123.jpg or 2345.bmp where the numbers will be the image's id's in the database. Don't repeat ids. These images will be deleted once added.

3) Create a blank file called 'recluster' (no quotes) and copy to IO dir. This will [re]create the bag-of-words dictionary.

4) To query the db, create a file called '123' (or whatever the id is) and put in IO. A file called 'matches' will appear (in the main folder)--this contains a list of id,match_strength pairs. All images in the db are returned. The first few are the ones appearing most similar. Talk to me about what the match strength actually means (a bit arbitrary).

5) To get feature matches (correspondences) between a pair of images (say 23 and 24) save a file called 'bb23,24' in the IO dir. a file called 'correspondences' will appear in the main folder. I will upgrade this to discard outliers at some stage...

Repeat 2,3,4 as needed. 'quit' will quit.

CONFIGURATION

The config file 'BoWSLAM/DefaultPatch.cfg' will be used by default. otherwise pass in the filename of a config file as a command line argument.

Config documentation will be coming soon, check out www.hilandtom.com/tombotterill/bow.