此流程仅用于 10x 平台单细胞准转录组、空间组和免疫组学数据分析。

## 流程图示例
![流程图](./config/dag.svg "流程图示例")
## 流程环境
``conda activate /public/home/weiyifan/miniforge3/envs/sk8``
## 流程部署
``snakedeploy deploy-workflow https://github.com/joyling93/sc10x_wf . --tag v1.1.4``
## 配置信息
config.yaml(config_space.yaml、config_multi为测试用例);
samples.yaml
## 流程运行
``snakemake -c30 --use-conda --cache``