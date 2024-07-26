# Chou-Fasman Method
Chou Fasman Algorithm to Predict Secondary Structural Elements (SSE)
# Backgroun Information
The Chou–Fasman method is an empirical technique for the prediction of secondary structures in proteins, originally developed in the 1970s by Peter Y. Chou and Gerald D. Fasman.[1][2][3] The method is based on analyses of the relative frequencies of each amino acid in alpha helices, beta sheets, and turns based on known protein structures solved with X-ray crystallography. From these frequencies a set of probability parameters were derived for the appearance of each amino acid in each secondary structure type, and these parameters are used to predict the probability that a given sequence of amino acids would form a helix, a beta strand, or a turn in a protein. The method is at most about 50–60% accurate in identifying correct secondary structures,[4] which is significantly less accurate than the modern machine learning–based techniques.[5]
# Comments and Steps
 I created my own HMM algorithm here to predict emission and transition values. I used the training dataset here but it’s incredibly long and my editor and computer can not always enough to give healthy outputs so I use  one thousands of the data. Especially, I’ve just used a little bit smaller data to train here to process my algorithm clearly sir. Also I researched some libraries like hmmlearn in python and look at the function’s features to get to know and understand more. In the propensity table  I’ve accepted some letters as if they would represent those that begin with same initial letter. (Like Glycine’s value for G). Here is the table.![image](https://github.com/user-attachments/assets/cc6428cc-ba18-4fa4-a891-7be237112c3a)
