# shared-analysis-tools
General-purpose code expected to be useful for multiple group members.

Install in your active environment with `pip install .`.

## IMPORTANT: `config.ini`
The user MUST put a config.ini in their home directory for certain methods that save data in this package to work. The contents should be
```
[paths]
figures = /path/to/which/the/code/should/save/figures
```

## Using `refactor` script
Use `refactor` in the terminal by navigating to your project folder and then executing 

~~~
UCI-tools.refactor -o old_text -n new_text
~~~
This will search all .py and .ipynb files in the project folder for `old_text`. `refactor` will give the user a preview of each instance where it finds `old_text` and ask for confirmation before replacing it with `new_text`.
