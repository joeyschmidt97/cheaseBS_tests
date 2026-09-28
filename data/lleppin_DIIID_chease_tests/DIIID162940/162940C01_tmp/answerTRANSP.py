#!/usr/bin/env python
import pexpect, os, sys, re


def b2s(obj):
    if isinstance(obj, bytes):
        for error_handling in ['strict', 'replace', 'ignore']:
            try:
                return obj.decode("utf-8", errors=error_handling)
            except Exception:
                pass
        raise RuntimeError('A bytes object was passed to b2s, but not handled')
    import numpy as np

    if isinstance(obj, np.ndarray) and (
        obj.dtype.name.startswith('bytes')
        or (obj.dtype.name.startswith('object') and np.all(map(lambda x: isinstance(x, (bytes, str)), obj.flat)))
    ):
        return np.reshape(np.array(list(map(b2s, obj.flat))), obj.shape)
    else:
        return obj


print('SHOT: ' + sys.argv[1])
shot = sys.argv[1]
print('DEVICE: ' + sys.argv[2])
device = sys.argv[2]
print('TRANSP version: ' + sys.argv[3])
version = ''
if sys.argv[3] != 'pshare':
    version = sys.argv[3]
print('MPI processors: ' + sys.argv[4])
trmpi = sys.argv[4]
print('NUBEAM MPI processors: ' + sys.argv[5])
ptrmpi = sys.argv[5]
print('TORIC MPI processors: ' + sys.argv[6])
toricmpi = sys.argv[6]
print('GENRAY MPI processors: ' + sys.argv[7])
genraympi = sys.argv[7]
print()

# =====================
# TRDAT
# =====================
if device == 'KSTR':
    pass
else:
    cases = ['\n', 'TRDAT TOP NODE >> ENTER OPTION CODE', pexpect.EOF]

    if os.system('which trdat') == 0:
        print('')
        print('#' * 20)
        print('# trdat')
        print('#' * 20)
        if not '#!/bin/bash' in open('env.sh', 'r').read():
            p = pexpect.spawn(
                os.environ.get('SHELL', '/bin/tcsh'), ['-c', 'source env.sh;setenv EDITOR cat;trdat ' + device + ' ' + shot], timeout=30
            )
        else:
            p = pexpect.spawn('/bin/bash', ['-c', 'source env.sh;export EDITOR=cat;trdat ' + device + ' ' + shot], timeout=30)
        while True:
            p.expect(cases)
            if p.after is pexpect.EOF:
                sys.stdout.write(b2s(p.before))
                break
            sys.stdout.write(b2s(p.before) + b2s(p.after))
            if 'TRDAT TOP NODE >> ENTER OPTION CODE' in b2s(p.after):
                p.sendline('Q')
            elif 'USER ACKNOWLEDGE - HIT ANY KEY' in b2s(p.after):
                break
            elif b2s(p.after) == '\n':
                pass

# =====================
# TR_START
# =====================
cases = [
    'tokyy_id: enter tokamak id:',
    'tokYY_id: verify run tok\.yy = \s+\.\s \(Y/N\):',
    'tr_start: verify your email address = .+@.+\..+ \(Y/N\):',
    'tr_start: Do you want MDSplus Output\? Y/N',
    'tr_start: Do you want to send Ufiles as tar file\? Y/N',
    'Press ENTER or type command to continue',
    '\[Ufiles USER ACKNOWLEDGE - HIT ANY KEY\]',
    'start: you are asking for NBI MPI',
    'start: you are asking for NB MPI',
    'start: you are asking for PTR MPI',
    'start: you are asking for TORIC MPI',
    'start: you are asking for GENRAY MPI',
    'Enter number of proccesses ',
    'btw 00 and 99',
    '\(Y/N\):',
    '     Verify Tree:',
    '\n',
    pexpect.EOF,
]

print('#' * 20)
print('# tr_start')
print('#' * 20)
if not '#!/bin/bash' in open('env.sh', 'r').read():
    p = pexpect.spawn(
        os.environ.get('SHELL', '/bin/tcsh'), ['-c', 'source env.sh;setenv EDITOR cat;tr_start ' + shot + ' ' + version], timeout=120
    )
else:
    p = pexpect.spawn('/bin/bash', ['-c', 'source env.sh;export EDITOR=cat;tr_start ' + shot + ' ' + version], timeout=120)

while True:
    p.expect(cases)
    if p.after is pexpect.EOF:
        sys.stdout.write(b2s(p.before))
        break
    sys.stdout.write(b2s(p.before) + b2s(p.after))
    if re.search(b2s(p.after), 'tokyy_id: enter tokamak id:'):
        p.sendline(device.upper())
    elif re.search(b2s(p.after), 'btw 00 and 99'):
        p.sendline('60')
    elif re.search(b2s(p.after), 'start: you are asking for NB MPI') or re.search(b2s(p.after), 'start: you are asking for NBI MPI'):
        p.sendline(str(max([4, int(trmpi)])))
    elif re.search(b2s(p.after), 'start: you are asking for PTR MPI'):
        p.sendline(ptrmpi)
    elif re.search(b2s(p.after), 'start: you are asking for TORIC MPI'):
        p.sendline(toricmpi)
    elif re.search(b2s(p.after), 'start: you are asking for GENRAY MPI'):
        p.sendline(genraympi)
    elif re.search(b2s(p.after), '     Verify Tree:'):
        p.sendline('transp_' + device.lower())
    elif b2s(p.after) == '\n':
        pass
    elif 'USER ACKNOWLEDGE - HIT ANY KEY' in b2s(p.after):
        break
    else:
        p.sendline('Y')
