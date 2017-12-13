#####################################################################################################
#   GIFT: Guided and Interpretable Factorization for Tensors - Applications to Human Cancer Analytics
#   
#   Authors: Sejoon Oh (ohhenrie@snu.ac.kr), Seoul National University
#            Jungwoo Lee (muon9401@gmail.com), Seoul National University
#            Lee Sael (sael@sunykorea.ac.kr), The State University of New York (SUNY) Korea
#   
#   Version : 1.0
#   Date: 2017-11-30
#   Main contact: Sejoon Oh
#
#   This software is free of charge under research purposes.
#   For commercial purposes, please contact the author.
######################################################################################################

1. Introduction

    GIFT is a Guided and Interpretable Factorization for Tensors. GIFT provides interpretable factor matrices by 
    encoding prior knowledge as a regularization term in its objective function.
    
    Please refer to the following website for the details of GIFT.
  	https://github.com/leesael/GIFT

2. Usage
	
	[Step 1] Install Armadillo, LAPACK, and BLAS libraries.

		GIFT requires Armadillo and OpenMP libraries.

		Armadillo library is available at the lib directory or http://arma.sourceforge.net/download.html.

		Notice that Armadillo needs LAPACK and BLAS libraries, and they are also available at the lib directory.

		OpenMP version 2.0 or higher is required for GIFT. (It is installed by default if you use gcc/g++ compiler)
 	
 	[Step 2] Compile and run GIFT.

		If you successfully install all libraries, "make stf, ptf, or gift" command will create the corresponding executable file.

		The executable file takes 6 (GIFT and Silenced-TF) or 5 (P-Tucker) arguments, which are the path of an input tensor file, path of a mask matrix (GIFT and Silenced-TF only), path of a directory for storing results, tensor order (N), lambda, and tensor ranks (N elements). The arguments MUST BE valid and in the above order.

		ex1) ./gift input.txt mask.txt result/ 3 0.01 10 20 30
		ex2) ./ptf input.txt result/ 3 0.01 10 20 30

		Note that an input tensor and mask matrices use base-1 indexing. Mask matrices should be given in a single file, where the first column indicates a mask matrix number and the rest of columns are indices in a mask matrix.
		
		If you put the command properly, GIFT and other methods will write all values of factor matrices and a core tensor in the result directory set by an argument. (PLEASE MAKE SURE THAT YOU HAVE A WRITE PERMISSION TO THE RESULT DIRECTORY!)

		ex) result/FACTOR1, result/CORETENSOR

3. Demo
	
	"make demo" command will run the demo program, which performs GIFT for a sample tensor with sample mask matrices.

	You can see factorization results in a 'sample' directory, while the intermediate process is presented on your screen.
	
	Note that you can check the sample tensor and mask matrices in the 'sample' directory and the descriptions are as follows.

	- sample.data: The sample tensor which is a 3-order tensor of size 4,555 (patient) * 3,994 (gene) * 5 (experiment type) with 613,355 nonzeros. 
	               Each row in sample.data corresponds to each nonzero. First three columns and the last column indicate a coordinate and a value of a nonzero, respectively. The indexing starts from 1. For example, a line ‘2 3 1 1.5’ means that the second patient’s third gene has value 1.5 at the first experiment type. 

	- sample.mask: The mask matrix M^(2) which is a matrix with size 3994 (gene) * 50 (number of gene sets, or rank) with 6,789 nonzeros. 
				   Each row in sample.mask corresponds to each intended entry. The first column indicates a mode, 2 for this sample mask. The second and the third column indicates a coordinate of an intended entry. The indexing starts from 1. For example, a line ‘2 3 40’ indicates an intended entry in the second factor matrix which is third in gene and 40th in rank. 



