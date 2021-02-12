/* the name of the folders need to be sample_1-2
 *  		first number is the number corresponding to the condition tested 
 *  		second number is the number corresponding to the repetition
 *  		the numbers are 0, 1, 2, 3, ...., 10, 11 (and not 00, 01, 02, 03)
 *  In the folder sample_*-*, there is a "stack" folder, we need the first frame named frame-0000.png
 */

// path and parameters to modify
path="D:/simulations/param_study/";
nb_sample=2;
nb_repetition=2;


//code it self
i_final=nb_sample-1;
k_final=nb_repetition-1;
File.makeDirectory(path+"stacks"); //create a folder to store .tif sequences 

for (i=0;i<=i_final;i++){
	for (k=0; k<=k_final;k++){
		run("Image Sequence...", "open="+path+"sample_"+i+"-"+k+"/stack/frame-0000.png sort use");
		saveAs("Tiff", path+"stacks/stack_"+i+"-"+k+".tif");	
		close();
	}
}
