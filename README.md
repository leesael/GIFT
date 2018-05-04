# GIFT

Overview
---------------

**Motivation**: Given multi-platform genome data with prior knowledge of functional gene sets, how can
we extract interpretable latent relationships between patients and genes? More specifically, how can we
devise a tensor factorization method which produces an interpretable gene factor matrix based on functional
gene set information while maintaining the decomposition quality and speed?

**Method**: We propose **GIFT**, a **G**uided and **I**nterpretable **F**actorization for **T**ensors. **GIFT** provides interpretable
factor matrices by encoding prior knowledge as a regularization term in its objective function.

**Results**: We apply GIFT to the PanCan12 dataset (TCGA multi-platform genome data) and compare the
performance with P-Tucker, our baseline method without prior knowledge constraint, and Silenced-TF, our
naive interpretable method. Results show that GIFT produces interpretable factorizations with high scalability
and accuracy. Furthermore, we demonstrate how results of GIFT can be used to reveal significant
relations between (cancer, gene sets, genes) and validate the findings based on literature evidence.

![overview_img](/img/overall.png)


Paper
---------------

**GIFT: Guided and Interpretable Factorization for Tensors with an Application to Large-Scale Multi-platform Cancer Analysis**  
[Sejoon Oh*](https://www.sejoonoh.com/), [Jungwoo Lee*](https://datalab.snu.ac.kr/~ljw9111/), and [Lee Sael](http://www3.cs.stonybrook.edu/~sael/) (* These authors contributed equally to this work.)   
[[Paper](/paper/GIFT.pdf)] [[Supplementary Material](/paper/Supplementary.pdf)]

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