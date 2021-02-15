# Residential Battery Management Tool (RBMT) [Behind the Meter BESS Management Strategies]

![Picture1](https://user-images.githubusercontent.com/69669859/97017890-5ef5e780-1546-11eb-9ec9-2bbfa502331a.jpg)


This repository contains a tool of three different power management strategies for the domestic residential batteries - Created October 2020 .

This tool can be used to generate the power dispatch of residential batteries (with any specifications) to minimize the household's electricity bill for any time series data (single day to multiple years) with any temporal resolution. 

The outputs of the RBMT are:
1.	The net household demand with and without the battery.
2.	Electricity bill with and without the battery or the PV.
3.	Battery power dispatch.
4.	Battery state of charge.
5.	Battery degradation. 
6.	Household’s voltage.
7.	Household’s losses. 
8.  PV self-consumption  
9.  Self-sufficiency 
10. Exported energy 
11. Curtailed energy  


Please check the RBMTGuide.pdf file for more details and guidance on how to use the code. 

This tool was validated and detailed in the following accepted paper, please acknowledge any contributions of the RBMT by citing: 

[1]. A. A. R. Mohamed, R. J. Best, X. Liu, and D. J. Morrow, ‘Domestic Battery Power Management Strategies to Maximize the Profitability and Support the Network’, IEEE PES General Meeting, pp.1-5, 2021


This code has been developed by [Ahmed A.Raouf Mohamed](https://pure.qub.ac.uk/en/persons/ahmed-mohamed) - EPIC Research Cluster, School of Electronics, Electrical Engineering and Computer Science at Queen's University in Belfast, UK. This work is part of [SPIRE 2 Project](https://www.ulster.ac.uk/spire2/the-project). 

For any inquiry: amohamed06@qub.ac.uk / AARaoufM@gmail.com 
[![twitter2](https://user-images.githubusercontent.com/69669859/97111234-a068cd00-16d5-11eb-9559-ff4b8946c0d8.png)](https://twitter.com/RA2OOOF)

v1.0 First release (10/2020).

v1.1 Added a selectable starting point of the BESS SOC (11/2020).

v1.2 Two modifications have been added to the first algorithm (CRBA): a) An option to charge part of the battery capacity overnight to optimize the time of use tariff; b) Adjusted the discharging to start after the end of low tariff period (11/2020).

v1.3 Standing charge has been added as per UK tariff structures (12/2020).

v1.4 Input data to be entered in a csv file (01/2021).

v1.5 Important version with the following features (02/2021): 
- The export power limit has been added to be entered as a utility input. 
- The PV self-consumption and the household self-sufficiency with and without the battery are calculated and printed.
- The exported energy with and without the battery is calculated and printed.
- The curtailed energy with and without the battery is calculated and printed.

Copyright @ 2020 
