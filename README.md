# 


<a name="readme-top"></a>



<h3 align="center">ngsAMOVA</h3>

  <p align="center">
    Genotype likelihood framework for Analysis of Molecular Variance with NGS data
    <br />
    <a href="https://github.com/github_username/repo_name"><strong>Quickstart»</strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">Installation</a>
    .
    <a href="https://github.com/github_username/repo_name">Tutorial</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#How to cite">How to cite</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## ngsAMOVA: [A]nalysis of [MO]lecular [VA]riance with NGS data



<!-- GETTING STARTED -->


 

### Installation

```
git clone https://github.com/isinaltinkaya/ngsAMOVA.git; cd ngsAMOVA; make
```


#### Prerequisites

Install dependencies
``` 
sudo apt-get install libz-dev liblzma-dev libbz2-dev libcurl4-openssl-dev g++ make
```
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage


_For more examples, please refer to the [Documentation](https://example.com)_




### Input files

#### Metadata file

If `-f/--formula` is not defined, ngsAMOVA will assume input metadata file is ordered and all of its hierarchical levels will be used.
```
./ngsAMOVA --formula "Individuals ~ Regions/Populations"
```


```
Individuals,Regions,Populations
ind1,reg1,pop1
ind2,reg1,pop2
ind3,reg2,pop3
ind4,reg2,pop3
ind5,reg2,pop4
```



<p align="right">(<a href="#readme-top">back to top</a>)</p>








<!-- How to cite -->
## How to cite

<p align="right">(<a href="#readme-top">back to top</a>)</p>

