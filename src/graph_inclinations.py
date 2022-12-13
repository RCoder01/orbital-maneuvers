import matplotlib.pyplot as plt

def graph(incs: list[float]):
    plt.hist(incs, bins=tuple(range(145)))
    plt.title('Debris in Low Earth Orbit by Inclination', fontsize=12)
    plt.xlabel('Inclination')
    plt.ylabel('Number of Cataloged Objects')
    plt.yscale('log')
    plt.gcf().set_size_inches(6, 3)
    plt.savefig('graphs/leo_inclinations.png', bbox_inches='tight', dpi=300)
    plt.clf()

def create_incs(input_filename: str):
    import json
    with open(input_filename, mode='r', encoding='UTF-8') as f:
        objects = json.load(f)
    incs = list(map(lambda o: o['INCLINATION'], objects))
    return incs

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        try:
            incs = create_incs(sys.argv[1])
        except FileNotFoundError:
            print(f'File with name {sys.argv[1]} not found')
            quit()
    else:
        try:
            with open('../data/incs.txt', mode='r', encoding='UTF-8') as f:
                incs = f.readlines()
        except FileNotFoundError:
            print('File with name ../data/incs.txt not found, '
            'try running this script with a debris file name as an argument\n'
            'Debris files can be created with get.py')
            quit()
        incs = list(map(lambda s: float(s.strip()), incs))
    graph(incs)
