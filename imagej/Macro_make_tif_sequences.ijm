/* the name of the folders need to be sample_1-2
 *  		first number is the number corresponding to the condition tested
 *  		second number is the number corresponding to the repetition
 *  		the numbers are 0, 1, 2, 3, ...., 10, 11 (and not 00, 01, 02, 03)
 *  In the folder sample_*-*, there is a "stack" folder, we need the first frame named frame-0000.png
 */

print("Starting up ImageJ analysis");
// path and parameters to modify
path="/data/u0063694/simulations/laning/sims_mathilde/parameter_study/";
max_n = 1;//Set this to -1 if you just want to process all directories that start with sample_

//code it self
File.makeDirectory(path+"stacks"); //create a folder to store .tif sequences

list = getFileList(path);
list = Array.sort(list);
n_sample = 0;//Counter of actual sample directories

for (i = 0; i < list.length; i++) {
   dirname = list[i];
   fullpath = path+File.separator+dirname;//Convert the directory name to a full path

   if( File.isDirectory( fullpath  ) ){//Only interested in directories
      if(  startsWith( dirname, "sample_" ) ){//Only interested in directories starting with sample_
         index_string = dirname.substring(7,dirname.length()-1);
         print("Running in directory "+dirname);
         run("Image Sequence...", "open="+fullpath+"/stack/frame-0000.png sort use");
         saveAs("Tiff", path+"/stacks/stack_"+index_string+".tif");
         n_sample += 1;//Increment the counter
         if( max_n > 0  ){
            if( n_sample >= max_n ){
               exit( "Stopped before all folders are processed" );
            }
         }
      }
   }
}
