# D I S C L A I M E R

You are using a BETA version of the program metric_element_create.f, which is currently under development by D. A. Spong of the Fusion Energy Division, Oak Ridge National Laboratory. Please report any problems or comments to him. As a BETA version, this program is subject to periodic change and improvement without notice.

7/29/2010: A test was included to see if all the surfaces to the edge were included in the xbooz_xform calculation. i.e., if(rmncbh(1,nsd) .eq. 0.)... This can also be seen from running ncdump on the boozmn file and comparing ns_b to comput_surfs. If the last surface is not included, then there can be problems with the outer surface in boozmn since averages between surfaces are made. This particularly shows up when xmetric is used for the outer surface metric calculation. In this case, the outer surface area can be too low by a factor of ~4 due to the outer R and z being down by 1/2 (from averaging over a 0 plus the correct value).
