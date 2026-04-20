import sys
import subprocess
import re

"""
Appends an element to a list and returns the list
"""
def rappend(list, appendix):
    list.append(appendix)
    return list

"""
Uses 'curl -I' to determine whether a file exists at a given URL.
"""
def check_file_curl(url):
    # Get the HTTP header
    cr = subprocess.run(rappend("curl -I".split(), url), capture_output = True)
    # curl failed to run
    if cr.returncode: return False
    success_str = b"HTTP/[0-9.]* 200"
    return bool(re.match(success_str, cr.stdout))

"""
Checks whether a corresponding URL exists. Uses curl for maximum portability.
"""
def check_file_http(url):
    # Check that curl exists
    if subprocess.run("which curl".split(), capture_output = True).returncode:
        raise "check_file_http: No valid command to fetch HTTP headers. Requires curl."
    return check_file_curl(url)
