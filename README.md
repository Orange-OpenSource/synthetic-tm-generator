# synthetic-tm-generator

Synthetic traffic matrix generator based on two publications from SIGCOMM CCR:
- A. Nucci, A. Sridharan, and N. Taft. [*The problem of synthetically generating IP traffic matrices*](https://dl.acm.org/citation.cfm?id=1070876);
- M. Roughan. [*Simplifying the synthesis of internet traffic matrices*](https://dl.acm.org/citation.cfm?id=1096551).

## Installation

### With Docker

Run the following two commands:
```
docker build -t tm-generator .
docker run -it tm-generator example_template.yml output_file.yml 0 1
```

### Without Docker

You will need to install a few Python dependencies:
```
sudo apt install -y python-tk
sudo pip install -r requirements.txt
```

If you want to use the GLPK solver, you'll also need to install it:
```
sudo apt install -y glpk-utils
```


## Examples

To generate attacks from all nodes to a single target:
```
./generate.py --max-node-class 1 template.yml output_file.yml 0 1
```

To generate attacks from a single node to 3 targets:
```
./generate.py --max-node-class 1 template.yml output_file.yml 1 3
```

To limit the network to only the largest nodes (first class of nodes only):
```
./generate.py --max-node-class 1 template.yml output_file.yml 0 1
```

To define a different maximum link capacity (defaults to 100 Gbps):
```
./generate.py --max-link-capacity 20000 template.yml output_file.yml 0 1
```

To define a different mean attack load (defaults to 10 Gbps):
```
./generate.py --max-link-capacity 500000 template.yml output_file.yml 1 1
```


## License

This project is licensed under [the Apache License v2.0](LICENSE).
