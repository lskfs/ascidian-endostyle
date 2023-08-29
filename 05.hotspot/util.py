
def get_cmap(cmap=None, N=None, category=False, theme='light'):

    import re
    import matplotlib as mpl

    if theme == 'light':
        facecolor = 'w'
    elif theme == 'dark':
        facecolor = 'k'

    is_hex_color = lambda x: re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', x)
    if isinstance(cmap, str) and is_hex_color(cmap):
        cmap = [facecolor, cmap]

    if isinstance(cmap, list):
        if category:
            if N is None:
                N = len(cmap)
            cmap = mpl.colors.ListedColormap(cmap[:N], N=N)
        else:
            cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', 
                    cmap, N=256)
    elif isinstance(cmap, str):
        if not category:
            N = 256
        else:
            assert N is not None
        cmap = mpl.cm.get_cmap(cmap, N)
        if category and isinstance(cmap, mpl.colors.LinearSegmentedColormap):
            cmap = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
            cmap = mpl.colors.ListedColormap(cmap, N=N)
    return cmap

