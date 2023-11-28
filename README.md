# shared-analysis-tools
General-purpose code expected to be useful for multiple group members.

Install in your active environment with `pip install .`.

## Using `refactor` script
Use `refactor` in the terminal by navigating to your project folder and then executing 

~~~
UCI-tools.refactor -o old_text -n new_text
~~~
This will search all .py and .ipynb files in the project folder for `old_text`. `refactor` will give the user a preview of each instance where it finds `old_text` and ask for confirmation before replacing it with `new_text`.
