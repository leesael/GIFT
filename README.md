# GIFT

Overview
---------------

**Motivation**: Given multi-platform genome data with prior knowledge of functional gene sets, how can we extract interpretable latent
relationships between patients and genes? More specifically, how can we devise a tensor factorization method which produces an
interpretable gene factor matrix based on gene set information while maintaining the decomposition quality and speed?

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

**GIFT: Guided and Interpretable Factorization for Tensors - An Application to Large-Scale Multi-platform Cancer Analysis**  
[Sejoon Oh*](https://www.sejoonoh.com/), [Jungwoo Lee*](https://datalab.snu.ac.kr/~ljw9111/), and [Lee Sael](http://www3.cs.stonybrook.edu/~sael/) (* These authors contributed equally to this work.)   
[[Paper](/paper/GIFT.pdf)]

Code
---------------
Refer to the code directory in our repository or download the following zip file.
[[GIFT-v1.0](https://datalab.snu.ac.kr/data/GIFT/GIFT.zip)]

Dataset
---------------
| Name | Structure | Size | Number of Nonzeros | Download |
| :------------ | :-----------: | :-------------: |------------: |:------------------: |
| PANCAN12 tensor     | Patient - Gene - Experiment Type | 4,555 &times; 14,351 &times; 5 | 180M | [DOWN](https://datalab.snu.ac.kr/data/GIFT/total.zip) |
| Mask matrix, **M**<sup>(2)</sup>	    | Gene - Gene set | 14,351 &times; 50 | 7K | [DOWN](https://datalab.snu.ac.kr/data/GIFT/mask.zip) |

The mask file contains information about unmasked entries (genes in gene set).

Tested Environment
---------------
We tested our proposed method **GIFT** in a Linux Ubuntu 16.04.3 LTS machine equipped with an Intel Xeon E5-2630 v4 2.2GHz CPU and 512GB RAM.