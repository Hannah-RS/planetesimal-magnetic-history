The output files `iie_sucess_info_merged.csv` and `pallasite_sucess_info_merged.csv` collate the key parameters for each run in Sanderson et. al. 2026 *EPSL* submitted. The meanings of columns not described under `run_info.csv` or `run_results.csv` in the [METADATA](https://github.com/Hannah-RS/planetesimal-magnetic-history/blob/main/METADATA.md) file are as below. The index *i* is replaced by 1-3 for the IIE irons (1 = Techado, 2 = Colomera, 3= Miles) and 1-5 for the Main Group pallasites (1=Marjalahti, 2=Brenham, 3=Springwater, 4=Imilac, 5=Esquel) denoting each meteorite in the order of remanence acquisition. Also note that the analysis flags f1, f2, and f3 have different meanings for the two meteorite groups. For more details see the post-processing scripts in the [iie_pallasite_processes](https://github.com/Hannah-RS/planetesimal-magnetic-history/tree/main/iie_pallasite_processes) folder

## For IIE irons and Main Group pallasites
+ d*i*_low - lower limit on meteorite formation depth [km]
+ d*i*_up - upper limit on meteorite formation depth [km]
+ c*i* - whether the core is solidifying at the time of remanence acquisition [boolean]

## For the IIE irons only
+ f1 - True if the dynamo is active at each meteorite's $^{40}$Ar/$^{39}$Ar age [boolean]
+ f2 - True if the location cooling through 623K at each meteorite's $^{40}$Ar/$^{39}$Ar age has the correct cooling rate [boolean]
+ f3 - True if the model and experimental normalised paleointensities agree [boolean]

If (f3 ==True) & (f2==True) then a model run is compatible with the IIE irons' thermal and magnetic history.

## For Main Group pallasites only
+ t*i*_low - lower limit on the time of remanence acquisition [Ma after CAI formation]
+ t*i*_up - upper limit on the time of remanence acquisition [Ma after CAI formation]
+ f1 - True if there is a gap in dynamo generation [boolean]
+ f2 - True if the positions with the required cooling rates at 925K also have the required cooling rates at 623K [boolean]
+ f3 - True if the dynamo is on/off as required for each meteorite's formation depth and time of remanence acquisition [boolean]
+ gm/f_dlow - lower limit on formation depth for Glorieta Mountain (gm) or Finmarken (f) [km]
+ gm/f_dup - upper limit on formation depth for Glorieta Mountain (gm) or Finmarken (f) [km]
+ gm/f_tlow - lower limit on time of remanence acquisition for Glorieta Mountain (gm) or Finmarken (f) [km]
+ gm/f_tup - upper limit on time of remanence acquisition for Glorieta Mountain (gm) or Finmarken (f) [km]
+ gm/f_mag - is the dynamo active at the time of remanence acquisition for Glorieta Mountain (gm) or Finmarken (f) [boolean]
+ gm/f_core_solid - is the core solidifying at the time of remanence acquisition for Glorieta Mountain (gm) or Finmarken (f) [boolean]

If f3 == True then a model run is compatible with the Main Group pallasites' thermal and magnetic history.