import os
import subprocess
import re

class Commands():
    content_length_title = b"Content-Length"
    content_length_lower = b"content-length"
    def __init__(self):
        self.query = None
        self.content_str = b""
        self.command = None

def touch_command(url, file_name):
    return ["touch", file_name]

def curl_command(url, file_name):
    return ["curl", "--output", file_name, url]

def curl_query(url):
    cr = subprocess.run(["curl", "-I", url], capture_output = True)
    return (cr.returncode, cr.stdout)

def wget_query(url):
    wr = subprocess.run(["wget", "--spider", "--server-response", url], capture_output = True)
    return (wr.returncode, wr.stderr)

def wget_command(url, file_name):
    return ["wget", url]

def dl_urls(urls_map, commands):
    """
    Downloads the SEACLIM hindcast forcing URLs, which are passed as an
    argument from the result of get_data_urls.base_urls. Checks for the
    presence and length of the file on both the server and the local file
    system.
    """

    http_check_pattern = b"HTTP/[0-9.]+ 200"
    file_length_pattern = b"[0-9]+"
    content_length_pattern = commands.content_str+b": "+file_length_pattern
    
    for model in urls_map:
        if not os.path.exists(model):
            os.mkdir(model)
        os.chdir(model)
        for expt in urls_map[model]:
            if not os.path.exists(expt):
                os.mkdir(expt)
            os.chdir(expt)
            for var in urls_map[model][expt]:
                for url in urls_map[model][expt][var]:
                    file_name = os.path.basename(url)
                    # Check the file on the server
                    # Get the HTTP header
                    (return_code, response) = commands.query(url)

                    server_existence = not return_code and bool(re.search(http_check_pattern, response))
                    if not server_existence:
                        print(f"URL not found for model {model}, experiment {expt}, variable {var}: {url}")
                        # If one URL is not found, then it is likely others of
                        # the same variable in the current experiment for the
                        # current model will also not be. Break the loop over
                        # times to avoid too many failure messages
                        break

                    server_length = int(re.search(file_length_pattern, re.search(content_length_pattern, response).group()).group())
                    local_existence = os.path.exists(file_name)
                    
                    len_equal = (os.path.getsize(os.path.basename(url)) == server_length) if local_existence else False
                    if not local_existence or not len_equal:
                        subprocess.run(commands.command(url, file_name))
            os.chdir("..")
        os.chdir("..")

if __name__ == "__main__":
    from get_data_urls import base_urls
    import argparse

    parser = argparse.ArgumentParser(
        prog = "seaclim_downloader",
        description = "Download the SEACLIM hindcast forcing files."
        )
    parser.add_argument("-d", "--dry-run", action="store_true")
    parser.add_argument("-c", "--curl", action="store_true")
    parser.add_argument("-w", "--wget", action="store_true")

    args = parser.parse_args()

    real_command = wget_command if args.wget else curl_command
    commands = Commands()
    if args.wget:
        real_command = wget_command
        commands.query = wget_query
        commands.content_str = Commands.content_length_title
    else:
        real_command = curl_command
        commands.query = curl_query
        commands.content_str = Commands.content_length_lower

    commands.command = touch_command if args.dry_run else real_command
        
    dl_urls(base_urls(), commands)
