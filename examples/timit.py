import subprocess

def std_call():
    cmd = ['../build/echellesimulator', '-s', 'MaroonX', '-t' ,'500' ,'-o','tmp.fit', '-k','0' ,'--fiber' ,'1']

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    td = ''
    for line in p.stdout:
        # print(line)
        if 'Duration Tracing: ' in str(line):
            print(str(line))
            td = float(str(line)[22:28])
    p.wait()

    print(p.returncode)
    return td

if __name__ =='__main__':
    times = []
    for i in range(5):
        td = std_call()
        times.append(td)
    print(min(times))