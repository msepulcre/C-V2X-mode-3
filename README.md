# Analytical performance model of LTE-V2X Mode 3 scheduling based on adaptive spatial reuse of radio resources
This code implements in Matlab the analytical models of the communication performance of C-V2X or LTE-V Mode 3 described in the following paper:

Paper title: *LTE-V2X Mode 3 scheduling based on adaptive spatial reuse of radio resources*  
Paper authors: *Daniel Sempere-García, Miguel Sepulcre and Javier Gozalvez*  

In order to comply with our sponsor guidelines, we would appreciate if any publication using this code references the above-mentioned publication.

This model quantifies the PDR (Packet Delivery Ratio) that could be achieved with our proposed scheme as a function of the distance between the transmitter and the receiver. In order to model the PDR, the following four mutually exclusive errors present in C-V2X are quantified:  
  1)	Errors due to half-duplex transmissions (HD)  
  2)	Errors due to a received signal power below the sensing power threshold (SEN)  
  3)	Errors due to propagation effects (PRO)  
  4)	Errors due to packet collisions (COL)  
  
CV2XMode3.m is the main script you have to run to get the PDR curve as a function of the distance for a given set of parameters, and the probability of each of the four transmission errors.

Input parameters (shall be set manually in the beginning of the code):  
* alpha: traffic density in veh/km. Values tested: 120.  
* lambda: pNumber of packets transmitted per second per vehicle. Values tested: 10, 20 and 50.  
* Psen: sensing threshold (dBm). Values tested: -90.5  
* Pt: transmission power in dBm. Values tested: 23.  
* S: number of sub-channels per sub-frame. Values tested: 2 and 4.  
* B: packet size in bytes. Values tested: 190.  

Output metrics:
* PDR: Packet Delivery Ratio for different Tx-Rx distances    
* deltaHD: probability of packet loss due to half-duplex transmissions for different Tx-Rx distances  
* deltaSEN: probability of packet loss due to a received signal power below the sensing power threshold for different Tx-Rx distances  
* deltaPRO: probability of packet loss due to propagation effects for different Tx-Rx distances  
* deltaCOL: probability of packet loss due to packet collisions for different Tx-Rx distances  
