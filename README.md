# GIFT

Overview
---------------

**Motivation**: Given cancer genome data with auxiliary gene set information, how can we extract significant
relations between cancers, gene sets, and genes? How can we devise a tensor factorization method which
produces interpretable factor matrices while maintaining the decomposition quality and speed?

**Method**: We propose GIFT, a **G**uided and **I**nterpretable **F**actorization for **T**ensors. GIFT provides interpretable
factor matrices by encoding prior knowledge as a regularization term in its objective function.

**Results**: Experiment results demonstrate that GIFT produces interpretable factorizations with high scalability
and accuracy, while other methods lack interpretability. We apply GIFT to the PANCAN12 dataset,
and GIFT reveals significant relations between cancers, gene sets, and genes, such as influential gene
sets for specific cancer (e.g., interferon-gamma response gene set for ovarian cancer) or relations between
cancers and genes (e.g., BRCA cancer $ APOA1 gene and OV, UCEC cancers $ BST2 gene).

![overview_img](/img/overview.pdf)


Paper
---------------

**GIFT: Guided and Interpretable Factorization for Tensors - Applications to Human Cancer Analytics**  
[Sejoon Oh*](https://www.sejoonoh.com/), [Jungwoo Lee*](https://datalab.snu.ac.kr/~ljw9111/), and [Lee Sael](http://www3.cs.stonybrook.edu/~sael/) (* These authors contributed equally to this work.)   
[[PDF](/paper/GIFT.pdf)], [[Supplementary Material](/paper/supple.pdf)]

Code
---------------
Refer to the code directory in our repository or download the following zip file.
[[GIFT-v1.0](/code/GIFT1.0.zip)]


Data
---------------
| Name | Structure | Size | Number of Entries | Download |
| :------------ | :-----------: | :-------------: |------------: |:------------------: |
| PanCan12     | Patient - Gene - Platform | 4,555 &times; 14,351 &times; 5 | 183,211,020 | [DOWN](https://datalab.snu.ac.kr/data/SNeCT/pancan12_tensor.tar.gz) |
| Pathway    | Gene - Gene | 14,351 &times; 14,351 | 665,429 | [DOWN](https://datalab.snu.ac.kr/data/SNeCT/pathway_network.tar.gz) |
