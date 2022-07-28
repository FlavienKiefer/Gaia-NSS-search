# Gaia-NSS-search
 Searching binary companions or planet in Gaia DR3 
   
You have to run   
   
   main_GaiaNSS.py    yourfile.extension   
   
Here yourfile.extension is e.g. "gaianss_test.txt" (everything beyond "." is considered extension).  
This file **must contain** at least one column called 'name' and be in the rdb-format (i.e. delimiter='\t', first row=2), like :  
  
name  
\-\-\-\-  
Beta Pictoris  
HD 172555  
HR 10  
HIP 65426  
  
**Dependences to install :**  
numpy, astroquery, astropy, argparse  
   
In case of troubles (and there will be for sure), please contact me at flavien.kiefer@obspm.fr   
Also if you think of improvements ;-) 

FAQ: If you have a certificate error, you might have to run beforehand "export SSL_CERT_FILE=/usr/local/lib/python3.9/site-packages/pip/\_vendor/certifi/cacert.pem" or similar to inform astroquery where to find the certificate "cacert.pem". If you don't know, use command "locate cacert.pem"
