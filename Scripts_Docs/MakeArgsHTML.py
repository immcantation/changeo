"""
Automatically generates html of presto command line usage
"""

import subprocess
from collections import OrderedDict

presto_path = '/home/jason/workspace/igpipeline/changeo_v0.4'
tools = OrderedDict([('AnalyzeAa', []),
                     ('CreateGermlines', []),
                     ('DefineClones', ['bygroup', 'hclust']),
                     ('MakeDb', ['imgt', 'igblast']),
                     ('ParseDb', ['delete', 'seq','clip']),
                     ('SplitDb', ['count', 'group', 'sort'])])

with open('output/arguments.html', 'w') as doc_handle, \
     open('output/navigation.html', 'w') as nav_handle: 
    # Start navigation list
    nav_handle.write('<ul>\n')
    for k, v in tools.iteritems():
        nav_handle.write('<li><a href="commandline.php#%s">%s</a></li>\n' % (k, k))
        #doc_handle.write('<div id=%s></div><h2>%s</h2>\n' % (k, k))
        doc_handle.write('<a class="anchor" id=%s></a><h2>%s</h2>\n' % (k, k))
        cmd = '/'.join([presto_path, k + '.py'])
        main_msg = subprocess.check_output([cmd, '--help'])
        doc_handle.write('<pre>\n%s\n</pre>\n' % main_msg)
        for c in v:
            doc_handle.write('<h3>%s</h3>\n' % c)
            sub_msg = subprocess.check_output([cmd, c, '--help'])
            doc_handle.write('<pre>\n%s</pre>\n' % sub_msg)
    # End navigation list
    nav_handle.write('</ul>\n')

print 'Wrote %s' % nav_handle.name
print 'Wrote %s' % doc_handle.name
