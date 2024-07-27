# Chou-Fasman Method
Chou Fasman Algorithm to Predict Secondary Structural Elements (SSE)
# Background Information
The Chou–Fasman method is an empirical technique for the prediction of secondary structures in proteins, originally developed in the 1970s by Peter Y. Chou and Gerald D. Fasman.[1][2][3] The method is based on analyses of the relative frequencies of each amino acid in alpha helices, beta sheets, and turns based on known protein structures solved with X-ray crystallography. From these frequencies a set of probability parameters were derived for the appearance of each amino acid in each secondary structure type, and these parameters are used to predict the probability that a given sequence of amino acids would form a helix, a beta strand, or a turn in a protein. The method is at most about 50–60% accurate in identifying correct secondary structures,[4] which is significantly less accurate than the modern machine learning–based techniques.[5]
# Comments and Steps
 I created my own HMM algorithm here to predict emission and transition values. I used the training dataset here but it’s incredibly long and my editor and computer can not always enough to give healthy outputs so I use  one thousands of the data. Especially, I’ve just used a little bit smaller data to train here to process my algorithm clearly sir. Also I researched some libraries like hmmlearn in python and look at the function’s features to get to know and understand more. In the propensity table  I’ve accepted some letters as if they would represent those that begin with same initial letter. (Like Glycine’s value for G). Here is the table.![image](https://github.com/user-attachments/assets/cc6428cc-ba18-4fa4-a891-7be237112c3a)
 And by this part in my code, I listed the numbers of structures according to the intervals in SSE table. ![image](https://github.com/user-attachments/assets/06ff26eb-b804-427a-bee0-c3df5cb85935)
 ![image](https://github.com/user-attachments/assets/110d9e65-4050-4ae7-9ee9-fa348f84f6c1)
 I created these functions the process the Chou-Fasman method clearly. 
And constructing the Chou-Fasman method here by sliding window 1 in every loop. It calculates the probabilites of alpha-helix and beta sheet situations and looking for overlaping after that. Also is_turn function looking for turns
![image](https://github.com/user-attachments/assets/c096bcac-c8a8-486f-9159-61469150342e) ![image](https://github.com/user-attachments/assets/ad8cfd5d-3ed9-47fa-b4cd-5198f7940c64)
And finally it compares the SSE values with the values that my algorithm’s proccesed ones in confusion matrix and calculating the precision, recall, accuracy and F1-score metrices for each SS element here.![image](https://github.com/user-attachments/assets/981edaf8-4158-4282-ba50-06e47b3883c8)
Here is the filled confusion matrix as my output.	 ![image](https://github.com/user-attachments/assets/304c88c2-01d2-46ef-9ce2-e74cab26143d)
And my calculated values for each SS element individually.![image](https://github.com/user-attachments/assets/d7ca9285-0f92-4dd2-86ce-0cdb645d866c)
Based on this, I guess my algorithm can be further improved for alpha helix and beta sheet values.






