import sys, subprocess, pathlib

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]
threads = int(sys.argv[3])

cmd = ['mafft', '--auto', '--thread', str(threads), str(in_fasta)]

with open(out_fasta, 'w') as out:
    try:
        subprocess.run(cmd, stdout=out, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        sys.exit(1)