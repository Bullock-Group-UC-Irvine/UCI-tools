def save_var_latex(key, value):
    dict_var = {}

    file_path = paths.paper +  "data.txt"

    try:
        with open(file_path, newline="") as file:
            reader = csv.reader(file)
            for row in reader:
                dict_var[row[0]] = row[1]
    except FileNotFoundError:
        pass

    dict_var[key] = value

    with open(file_path, "w") as f:
        for key in dict_var.keys():
            f.write(f"{key},{dict_var[key]}\n")

    return None

def save_prediction(string, y, dy):
    '''
    Save strings to data.txt to be used by paper.tex.

    Parameters
    ----------
    string: str 
	The identifier for the variable
    y: str 
        The value of the target variable
    dy: np.ndarray of str
	The uncertainty in the variable. The user *must* provide it as an
        np.ndarray, even if it's only one single symmetric uncertainty.

    Returns
    -------
    None
    ''' 

    save_var_latex(string, y)
    if len(dy)==1:
        save_var_latex('d'+string, dy[0])
        save_var_latex('d'+string+'_plus', 'n/a')
        save_var_latex('d'+string+'_minus', 'n/a')
    elif len(dy)==2:
        save_var_latex('d'+string, 'n/a')
        save_var_latex('d'+string+'_plus', dy[0])
        save_var_latex('d'+string+'_minus', dy[1])
    else:
        raise ValueError('Margin of error has a maximum of two elements.')
    return None
