#!/usr/bin/env python3
"""Group ThreadSanitizer reports by first project frame: tsan_summary.py <log>

The frame pattern expects the source tree to sit in a directory named "tgg".
"""
import sys,re,collections
txt=open(sys.argv[1]).read()
reps=txt.split('==================')
groups=collections.Counter(); ex={}
def userframe(block):
    # first frame in project code
    for l in block:
        m=re.search(r'#\d+ (.*?) (\S*/tgg/(pangenie/src|src|PopMergeVCF)/\S+?):(\d+)',l)
        if m:
            fn=re.sub(r'\(.*','',m.group(1)); fn=re.sub(r'std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >','string',fn)
            return '%s @ %s:%s'%(fn.strip()[:70], m.group(2).split('/tgg/')[1], m.group(4))
    for l in block:
        m=re.search(r'#\d+ (.*?) (/usr/include/jellyfish/\S+|\S+libgomp\S*|\S+libjellyfish\S*)',l)
        if m: return (re.sub(r'<.*','',m.group(1))[:50]+' @ '+m.group(2).split('/')[-1])
    return block[0].strip()[:80] if block else '?'
for r in reps:
    if 'WARNING: ThreadSanitizer' not in r: continue
    kind=re.search(r'WARNING: ThreadSanitizer: (.*?) \(',r).group(1)
    lines=r.split('\n'); blocks=[]; cur=None; heads=[]
    for l in lines:
        if re.match(r'^  (Write|Read|Previous|Atomic|Location|Mutex|Thread)',l) or (l.startswith('  ') and not l.startswith('    ') and l.strip()):
            cur=[]; blocks.append((l.strip(),cur))
        elif l.startswith('    #') and cur is not None: cur.append(l)
    acc=[(h,b) for h,b in blocks if re.match(r'(Write|Read|Previous|Atomic)',h)]
    if len(acc)<2: key=(kind,'?','?')
    else:
        a=re.sub(r' at 0x\w+| by .*','',acc[0][0]); b=re.sub(r' at 0x\w+| by .*','',acc[1][0])
        ta='main' if 'main thread' in acc[0][0] else 'worker'; tb='main' if 'main thread' in acc[1][0] else 'worker'
        key=(kind, '%s[%s] %s'%(a,ta,userframe(acc[0][1])), '%s[%s] %s'%(b,tb,userframe(acc[1][1])))
    groups[key]+=1
    if key not in ex: ex[key]=r
for k,c in groups.most_common():
    print(c,'|',k[0]); print('     A:',k[1]); print('     B:',k[2])
print('distinct groups',len(groups),'reports',sum(groups.values()))
