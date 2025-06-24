
source setup_env.sh

if [ "$#" -eq 0 ]; then
    condorDir=$(ls ${MBD_CENTRALITY_CALIB_PATH}/macros/condorDir/ -Artl | tail -1 | cut -d ":" -f 2 | cut -d " " -f 2)
    echo "Using the most recent job."
else

    condorDir=${1}
    echo "Using ${condorDir}."
fi
    
dir=${MBD_CENTRALITY_CALIB_PATH}/macros/condorDir/${condorDir}
nfiles=$(ls ${dir}/*.out | wc -l)

echo "runnumber,events,ZDC,ZDCone,mb,vtxcut,vtx,vtxsig,fillpattern,mean,mean_bal,scale,cent_chi2,cent_chi2_bal,ns,chi2,mu,k,trigeff,trigeff_err,npart,ecc2,ecc3,b,chi2_bal,mu_bal,k_bal,trigeff_bal,trigeff_err_bal,npart_bal,ecc2_bal,ecc3_bal,b_bal" > table.csv

for i in $(seq 0 ${nfiles}) ; do cat ${dir}/condor_${i}.out | grep "Table Entry" -A 1 | tail -n 1 >> table.csv; done;
