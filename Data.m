%% Copy the data from Excel or CSV file in the following format
%Or read from CSV file directly as
%Profile=readmatrix('Data.csv');
% | Demand(kW) | PV (kW) | EV(kW) | 
Profile=[
0.239725856	0	0
0.203014791	0	0
0.151259353	0	0
0.222635868	0	0
0.094974736	0	0
0.208889396	0	0
0.102328095	0	0
0.10932551	0	0
0.112346395	0	0
0.01431772	0	0
0.203082893	0	0
0.106599387	0	0
0.042571134	0	0
0.096148069	0	0
0.185280796	0	0
0.075356284	0	0
0.328307454	0	0
0.392836262	0	0
0.128706365	0	0
0.101588065	0	0
0.12070557	0	0
0.047746893	0	0
0.24836522	0	0
0.027887888	0	0
0.127064709	0	0
0.05398466	0.044	0
0.129450631	0.049	0
0.261248852	0.069	0
0.121547678	0.082	0
0.07620347	0.127	0
0.1692957	0.112	0
0.228802076	0.136	0
0.263402223	0.099	0
0.114745157	0.072	0
0.215050296	0.055	0
0.007704548	0.105	0
0.130059765	0.123	0
0.064806505	0.1	0
0.12913373	0.11	0
0.103610741	0.178	0
0.005797319	0.267	0
0.21459194	0.321	0
0.196324538	0.31	0
0.243248913	0.648	0
0.127426549	1.011	0
0.227951077	1.101	0
0.238483309	1.232	0
0.127565614	2.002	0
0.099751706	2.226	0
0.216007866	2.088	0
0.008240428	2.772	0
0.178973076	2.66	0
0.119301873	2.681	0
0.025725052	2.803	0
1.137250733	2.771	0
0.566930924	2.879	0
1.130818343	2.94	0
2.289253589	2.997	0
0.668870252	3.316	0
1.267376158	3.409	0
0.108756414	1.475	0
1.359091809	0.662	0
3.041151777	0.522	0
1.007718409	3.517	0
0.41736017	3.312	0
0.721421421	3.289	0
0.058220342	3.432	0
0.114903715	1.981	0
0.150875943	2.192	0
0.232835034	0.978	0
0.038862066	0.961	0
0.0523029	1.806	0
0.566863179	3.52	0
0.492652628	0.744	0
1.545984193	0.986	0
0.24600973	2.035	0
0.282300674	0.72	0
0.025189595	3.484	0
0.11568452	3.36	0
0.339497644	3.172	0
0.233317836	3.093	0
0.724836336	3.122	0
0.414657023	2.904	0
0.318506641	2.807	0
0.420938446	2.746	0
0.261588237	2.644	0
0.410973318	2.57	0
0.48836902	2.48	0
0.071149902	2.404	0
0.182981078	1.774	0
0.124560224	2.159	0
0.243011205	2.032	0
0.293928571	1.935	0
0.216362087	1.833	0
0.200518905	1.756	0
0.150119008	1.612	0
0.780427702	1.457	0
0.449746523	1.352	0
1.037825775	1.137	0
0.294706121	1.033	0
0.364405119	0.897	0
0.407388759	0.744	0
0.026768846	0.61	0
0.623021487	0.448	0
0.200709668	0.259	0
0.278568482	0.141	0
0.299003162	0.133	0
0.272928356	0.121	0
0.124513194	0.11	0
1.691630091	0.098	3.696
1.167356714	0.094	3.696
0.303494043	0.087	3.701
1.074391246	0.081	3.701
0.25561471	0.074	3.701
0.363982587	0.066	3.706
1.471772187	0.058	3.706
1.269245226	0.051	3.706
0.022826157	0.038	3.706
0.580509087	0.029	3.711
0.557664756	0.014	3.727
0.047122725	0.002	3.666
0.338131515	0	2.149
0.73524576	0	1.118
0.064341734	0	0
0.49669298	0	0
0.626965286	0	0
0.734624554	0	0
0.459650653	0	0
0.101724793	0	0
0.131049719	0	0
0.38057484	0	0
0.149875441	0	0
0.244804549	0	0
0.033654721	0	0
0.05904073	0	0
0.005509421	0	0
0.183492939	0	0
0.14849764	0	0
0.111062367	0	0
0.185379566	0	0
0.324558068	0	0
0.34822852	0	0
0.03738072	0	0
0.100390761	0	0
0.034917023	0	0
0.172416425	0	0
0.116666552	0	0
0.066951311	0	0
0.224052955	0	0
0.032995734	0	0
0.001444615	0	0
0.197326698	0	0
0.125228687	0	0
0.250234556	0	0
0.158763869	0	0
0.293001575	0	0
0.086236634	0	0
0.001910005	0	0
0.303353361	0	0
0.100096764	0	0
0.044748156	0	0
0.17915508	0	0
0.047393479	0	0
0.209523973	0	0
0.067082548	0	0
0.130225403	0	0
0.127465732	0	0
0.066308865	0	0
0.078723027	0	0
0.33459084	0.034	0
0.356186133	0.057	0
0.07219932	0.069	0
0.12585733	0.067	0
0.12594335	0.124	0
0.10102216	0.169	0
0.138708325	0.396	0
0.084269515	0.64	0
0.099792256	0.715	0
0.114976884	0.662	0
0.12273086	0.566	0
0.02939973	0.749	0
0.140260294	1.124	0
0.167839976	1.061	0
0.161291604	1.202	0
0.374971925	1.217	0
0.246736472	1.241	0
0.103195571	1.476	0
0.160926655	1.799	0
0.059877774	1.973	0
1.009598829	2.064	0
0.620713555	2.142	0
1.056187615	2.325	0
2.269623648	2.413	0
1.817663491	2.26	0
3.378212861	2.322	0
0.378530578	2.376	0
0.454313048	1.869	0
0.031156375	2.245	0
1.072975376	2.751	0
0.229904589	2.357	0
1.140620036	2.561	0
0.029908552	2.713	0
0.113675597	3.103	0
0.19391585	3.048	0
0.154089536	2.547	0
0.073408225	2.668	0
0.096502239	3.333	0
0.156115576	3.312	0
0.274110639	1.87	0
0.082773785	3.506	0
0.264886793	3.503	0
0.015625704	1.03	0
0.313487503	2.583	0
0.047102635	3.145	0
0.128092404	0.992	0
0.148804961	3.425	0
0.04707193	3.41	0
0.095657454	3.52	0
0.181270616	2.968	0
0.801017014	1.861	0
0.140390585	3.298	0
1.340092401	3.386	0
0.141567136	2.981	0
0.513295545	0.817	0
0.29013732	0.965	0
0.599948224	2.459	0
0.047775456	2.098	0
0.28377632	2.149	0
0.746479893	1.882	0
1.197190816	2.029	0
0.553829291	1.926	0
0.207917646	1.784	0
0.030550067	1.83	0
0.463532287	1.697	0
0.032197169	1.493	0
0.660566623	1.485	0
0.319736208	1.429	0
2.212380469	1.117	0
2.358353763	1.016	0
1.733765768	1.109	0
0.989490945	1.026	0
2.220391627	0.952	0
1.596117428	0.92	0
1.124818297	0.876	0
0.963045183	0.811	0
1.50313652	0.746	0
1.262612991	0.596	0
0.641357855	0.566	0
0.256029154	0.499	0
0.925135395	0.486	0
0.067130344	0.452	0
0.020234261	0.417	0
0.442485677	0.364	0
0.133565839	0.226	0
0.436448484	0.178	0
0.396632013	0.152	0
0.101778829	0.123	0
0.244089158	0.081	0
0.380143572	0.054	0
0.279805176	0.039	0
0.096051251	0.029	0
0.156852147	0.029	0
0.086568558	0.009	0
0.134579295	0.002	1.557
0.146888614	0	3.706
0.099929591	0	3.706
0.077181796	0	3.706
0.045450852	0	3.706
0.124805271	0	3.706
0.153743877	0	3.691
0.243465383	0	3.609
0.091373512	0	3.686
0.259161105	0	3.681
0.20615983	0	3.711
0.324578429	0	3.711
0.171261741	0	3.711
0.257246556	0	3.717
0.154907273	0	3.722
0.154846171	0	3.717
0.237399103	0	3.742
0.042312154	0	3.752
0.314288743	0	3.737
0.248641999	0	3.788
0.209323688	0	3.798
0.190034313	0	1.093
0.304292437	0	0
1.182457408	0	0
0.686750155	0	0

];