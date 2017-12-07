# GIFT

Overview
---------------

**Motivation**: Given cancer genome data with auxiliary gene set information, how can we extract significant
relations between cancers, gene sets, and genes? How can we devise a tensor factorization method which
produces interpretable factor matrices while maintaining the decomposition quality and speed?

**Method**: We propose **GIFT**, a **G**uided and **I**nterpretable **F**actorization for **T**ensors. **GIFT** provides interpretable
factor matrices by encoding prior knowledge as a regularization term in its objective function.

**Results**: Experiment results demonstrate that **GIFT** produces interpretable factorizations with high scalability
and accuracy, while other methods lack interpretability. We apply **GIFT** to the PANCAN12 dataset,
and **GIFT** reveals significant relations between cancers, gene sets, and genes, such as influential gene
sets for specific cancer (e.g., interferon-gamma response gene set for ovarian cancer) or relations between
cancers and genes (e.g., BRCA cancer *<->* APOA1 gene and OV, UCEC cancers *<->* BST2 gene).

![overview_img](/img/overall.png)


Paper
---------------

**GIFT: Guided and Interpretable Factorization for Tensors - Applications to Human Cancer Analytics**  
[Sejoon Oh*](https://www.sejoonoh.com/), [Jungwoo Lee*](https://datalab.snu.ac.kr/~ljw9111/), and [Lee Sael](http://www3.cs.stonybrook.edu/~sael/) (* These authors contributed equally to this work.)   
[[Paper](/paper/GIFT.pdf)], [[Supplementary Material](/paper/supple.pdf)]

Code
---------------
Refer to the code directory in our repository or download the following zip file.
[[GIFT-v1.0](/code/GIFT1.0.zip)]

Dataset
---------------
| Name | Structure | Size | Number of Nonzeros | Download |
| :------------ | :-----------: | :-------------: |------------: |:------------------: |
| PANCAN12 tensor     | Patient - Gene - Experiment Type | 4,993 &times; 14,591 &times; 5 | 180M | [DOWN](https://datalab.snu.ac.kr/data/GIFT/total.zip) |
| Mask matrix, **M**<sup>(2)</sup>	    | Gene - Gene set | 14,591 &times; 50 | 7K | [DOWN](https://datalab.snu.ac.kr/GIFT/mask.zip) |

**M**<sup>(1)</sup> and **M**<sup>(3)</sup> are filled with zeros, and three mask matrices are concatenated together in a single mask file. The file contains information about intended entries (genes in gene set) which are to be initialized as zeros in **GIFT**.
