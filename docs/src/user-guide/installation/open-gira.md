# open-gira 

open-gira relies on many tools, some of which are complex to install. 
The installation process may take a while. 
Users familiar with Docker may wish to use the Docker install options to 
simplify the process.

## Notes:
- If you are trying to install on **Windows**, use Windows Subsystems for Linux (WSL).
- Linux users may find the version of osmium-tool provided by their package repository
is out of date.
The instructions for building osmium-tool from source on Linux should work, however.

## Tests
Once everything is installed, run the tests to ensure everything is linked up properly:

```
python -m pytest tests
```