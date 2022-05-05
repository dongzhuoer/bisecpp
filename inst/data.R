# clean
dir('data', full.names = T) |> file.remove()



# Process raw data
ideal.para11 <- c(K.i = 0.9, K.nat = 0.1, k.d = 0.02);        #protein (with IPTG)
ideal.para12 <- c(K.i = 0, K.nat = 0.1, k.d = 0.02);          #protein (without IPTG)
ideal.para2 <- c(k.in = 0.008, k.out = 0.045, k.p = 0.006);   #AI-2
usethis::use_data(ideal.para11, ideal.para12, ideal.para2, overwrite = TRUE)

extra11 <- c(mu = 0.0044, iOD600 = 0.5);                      #cell (with ACDB)
extra12 <- c(mu = 0.0056, iOD600 = 0.5);                      #cell (others)
extra2 <- c(IPTG.ACDB = TRUE, IPTG.K = TRUE, AI2.out.0 = 1);  #sundries
usethis::use_data(extra11, extra12, extra2, overwrite = TRUE);

space1 <- cbind(integer(length(ideal.para11)), ideal.para11 * 3);
space2 <- cbind(integer(length(ideal.para2)), ideal.para2 * 3);
usethis::use_data(space1, space2, overwrite = TRUE)


protein_ACDB <- bisecpp::f_protein(ideal.para11, extra11) # protein with ACDB gene and IPTG
protein_K <- bisecpp::f_protein(ideal.para11, extra12) # protein with IPTG
protein <- bisecpp::f_protein(ideal.para12, extra12) # protein with nothing
usethis::use_data(protein_ACDB, protein_K, protein, overwrite = TRUE);

data1 <- bisecpp::simulate_data(bisecpp::f_protein, ideal.para11, extra11, seq(60,240,60));
data2 <- bisecpp::simulate_data(bisecpp::f_AI2_out, ideal.para2, extra2, seq(0,270,30));
usethis::use_data(data1, data2, overwrite = TRUE);
