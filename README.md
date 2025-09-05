# UCI-tools
General-purpose code expected to be useful for multiple group members.

Install in your active environment with `pip install .`.

## `config.ini`
The user should put a `<environment>_config.ini` in their home directory for certain methods that save data in this package to work. The user must name it according to which conda environment they plan to use: `<environment>_config.ini`. If you use the `base` environment (which you shouldn't), you can name the config file `config.ini`. The contents should be
```
[paths]
figures = /path/to/which/the/code/should/save/figures
```
If you don't create this before the first time you import `UCI_tools`, the package will create one for you, and it will default to saving all outputs in `$HOME/<environment>_output/`.

## Using `refactor` script
Use `refactor` in the terminal by navigating to your project folder and then executing 

~~~
UCI-tools.refactor -o old_text -n new_text
~~~
This will search all .py and .ipynb files in the project folder for `old_text`. `refactor` will give the user a preview of each instance where it finds `old_text` and ask for confirmation before replacing it with `new_text`.
