# Vasptools
some toolkits for VASP
## vaspeqstress.sh
小侯飞氚用python2写了VASP定压优化脚本，我之前用python3重写成vaspeqstress.py，但是我不喜欢python加模块的模式，而且服务器上常常不好自己配置python，当然也可能是我想多了。总之我用bash重新写了vaspeqstress.sh，不需要任何多余的模块，并且这样也可以很方便的用于作业管理系统。
## 晶格变形脚本的区别
deform.sh用于向POSCAR文件施加变形，seriesdeform.sh用于施加系列变形，idealdeform用于计算应力应变曲线。注意，以上三个脚本施加的变形都是基于晶体坐标系施加变形，也就是晶格的a b c分别对应[1 0 0] [0 1 0] [0 0 1]。另外三个带后缀2的脚本功能与前者完全相同，但是变形是在笛卡尔坐标系上，跟晶格坐标系没有关系。
