#! /bin/bash

# KofamScan uses KEGG Orthologs to assign functions to predicted ORFs
# First, input file with ORFs predicted previously using prodigal
# -o : for output file
# -f : mapper to use functional profiling
# -E : e-value max accepted
# -k : official list of prokaryotic KOs from KEGG database
# -p : functional profiles from KEGG database
# --cpu : number of cpu used

exec_annotation AAP38-BRE-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/BRE-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-BRE-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/BRE-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-BRE-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/BRE-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-BRE-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/BRE-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-GAB-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/GAB-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-GAB-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/GAB-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-GAB-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/GAB-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-GAB-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/GAB-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-MAT-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MAT-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MAT-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MAT-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MAT-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MAT-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MAT-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MAT-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-MGA-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MGA-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MGA-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MGA-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MGA-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MGA-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MGA-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MGA-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-MOL-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MOL-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MOL-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MOL-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MOL-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MOL-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-MOL-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/MOL-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-PON-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PON-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-PON-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PON-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-PON-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PON-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-PON-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PON-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-PLA-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PLA-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-PLA-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PLA-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-PLA-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PLA-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-PLA-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/PLA-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-SAU-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/SAU-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-SAU-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/SAU-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-SAU-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/SAU-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-SAU-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/SAU-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8

exec_annotation AAP38-TON-A-S2_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/TON-A-S2_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-TON-A-S3_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/TON-A-S3_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-TON-A-S4_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/TON-A-S4_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8
exec_annotation AAP38-TON-A-S5_output_protein.faa -o /home/serveur-mhugoni/Bureau/2026_Flora/TON-A-S5_kofam.txt -f mapper -E 1e-15 -k ko_list_07_2025 -p profiles_07_2025 --cpu 8