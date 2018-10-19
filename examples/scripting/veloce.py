import subprocess

cmd = ['../../build/echellesimulator', '-s', 'VeloceRosso', '-t', '0.1', '-o', 'veloce_mdwarf.fit', '-k', '0',
       '--fiber', '10', '--phoenix', '3200,4.50,-0.0,-0.0,0']

p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
for line in p.stdout:
    print(line)
p.wait()

print(p.returncode)

for i in range(11, 13):
    cmd = cmd = ['../../build/echellesimulator', '-s', 'VeloceRosso', '-t', '0.1', '-o', 'veloce_mdwarf.fit', '-k', '1',
                 '--fiber', str(i), '--phoenix', '3200,4.50,-0.0,-0.0,0']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in p.stdout:
        print(line)
    p.wait()

    print(p.returncode)
