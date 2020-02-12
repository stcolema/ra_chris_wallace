#!/usr/env/bash

# Common variables across calls
CMD="/home/stephen/Documents/PhD/Year_1/Consensus_clustering/ra_chris_wallace/Gen_data/generate_data_pk.R";
N_CLUST="25 50 75 100 150";

# Call specific variables for near clusters
MEANS="1 1.5 3 4 5";
P=100;
FILENAME="near_clusters.csv";
PH_SAVENAME="near_clusters_ph.png";
COR_PH_SAVENAME="near_clusters_cor_ph.png";

Rscript "${CMD}" -p "${P}" --n_clust "${N_CLUST}" --means "${MEANS}" --file ${FILENAME} --ph_save_name ${PH_SAVENAME} --cor_ph_save_name ${COR_PH_SAVENAME}

P=250;
FILENAME="near_clusters_p_250.csv";
PH_SAVENAME="near_clusters_ph_p_250.png";
COR_PH_SAVENAME="near_clusters_cor_ph_p_250.png";

Rscript "${CMD}" -p "${P}" --n_clust "${N_CLUST}" --means "${MEANS}" --file ${FILENAME} --ph_save_name ${PH_SAVENAME} --cor_ph_save_name ${COR_PH_SAVENAME}

# Call specific variables for mediumly spaced clusters
MEANS="1 2 4 6 8";
P=100;
FILENAME="medium_clusters.csv";
PH_SAVENAME="medium_clusters_ph.png";
COR_PH_SAVENAME="medium_clusters_cor_ph.png";

Rscript "${CMD}" -p "${P}" --n_clust "${N_CLUST}" --means "${MEANS}" --file ${FILENAME} --ph_save_name ${PH_SAVENAME} --cor_ph_save_name ${COR_PH_SAVENAME}

P=250;
FILENAME="medium_clusters_p_250.csv";
PH_SAVENAME="medium_clusters_ph_p_250.png";
COR_PH_SAVENAME="medium_clusters_cor_ph_p_250.png";

Rscript "${CMD}" -p "${P}" --n_clust "${N_CLUST}" --means "${MEANS}" --file ${FILENAME} --ph_save_name ${PH_SAVENAME} --cor_ph_save_name ${COR_PH_SAVENAME}

# Call specific variables for distant clusters
MEANS="1 3 6 9 12";
P=100;
FILENAME="distant_clusters.csv";
PH_SAVENAME="distant_clusters_ph.png";
COR_PH_SAVENAME="distant_clusters_cor_ph.png";

Rscript "${CMD}" -p "${P}" --n_clust "${N_CLUST}" --means "${MEANS}" --file ${FILENAME} --ph_save_name ${PH_SAVENAME} --cor_ph_save_name ${COR_PH_SAVENAME};

P=250;
FILENAME="distant_clusters_p_250.csv";
PH_SAVENAME="distant_clusters_ph_p_250.png";
COR_PH_SAVENAME="distant_clusters_cor_ph_p_250.png";

Rscript "${CMD}" -p "${P}" --n_clust "${N_CLUST}" --means "${MEANS}" --file ${FILENAME} --ph_save_name ${PH_SAVENAME} --cor_ph_save_name ${COR_PH_SAVENAME};
