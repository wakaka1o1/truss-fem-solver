# 1D Truss FEM Solver

一个 C++20 的 1D 平面桁架有限元求解器，JSON 输入，控制台输出节点位移、杆单元应力和支座反力。

## 项目背景

本项目是 [Superpowers](https://github.com/anthropics/claude-code) + Claude Code 的完整工作流测试 Demo，覆盖了从设计到交付的全流程：

| 阶段 | Superpowers Skill | 说明 |
|---|---|---|
| 需求设计 | `brainstorming` | 交互式对话确定架构、输入格式、技术选型 |
| 实现计划 | `writing-plans` | 10 个 TDD 任务，完整代码，零占位符 |
| 环境隔离 | `using-git-worktrees` | Git worktree 隔离开发分支 |
| 开发执行 | `subagent-driven-development` | 独立子代理逐任务实现 + 规范审查 + 代码审查 |
| 完成交付 | `finishing-a-development-branch` | 合并到 master，清理工作区 |

## 技术栈

- **语言：** C++20
- **构建：** CMake
- **测试：** Google Test
- **JSON 解析：** nlohmann/json
- **求解器：** 高斯消元（部分选主元）

## 构建与运行

```bash
# 构建
cmake -B build && cmake --build build

# 运行测试
ctest --test-dir build

# 运行示例
./build/src/Debug/truss_solver tests/fixtures/two_bar_truss.json
./build/src/Debug/truss_solver tests/fixtures/three_bar_truss.json
```

## 输入格式

```json
{
  "nodes": [
    { "id": 0, "x": 0.0,  "y": 0.0 },
    { "id": 1, "x": 1.0,  "y": 0.0 },
    { "id": 2, "x": 0.5,  "y": 1.0 }
  ],
  "elements": [
    { "id": 0, "n1": 0, "n2": 2, "E": 210000000000.0, "A": 0.001 },
    { "id": 1, "n1": 1, "n2": 2, "E": 210000000000.0, "A": 0.001 }
  ],
  "boundary_conditions": [
    { "node": 0, "dof": "ux", "value": 0.0 },
    { "node": 0, "dof": "uy", "value": 0.0 },
    { "node": 1, "dof": "ux", "value": 0.0 },
    { "node": 1, "dof": "uy", "value": 0.0 }
  ],
  "loads": [
    { "node": 2, "dof": "uy", "value": -10000.0 }
  ]
}
```

## 输出示例

```
=== Nodal Displacements ===
Node 0: ux = 0.000000e+00, uy = 0.000000e+00
Node 1: ux = 0.000000e+00, uy = 0.000000e+00
Node 2: ux = 0.000000e+00, uy = -3.327482e-05

=== Element Stresses ===
Element 0: sigma = -5.590151e+06 Pa
Element 1: sigma = -5.590151e+06 Pa

=== Support Reactions ===
Node 0: Fx = 2.500000e+03
Node 0: Fy = 5.000000e+03
Node 1: Fx = -2.500000e+03
Node 1: Fy = 5.000000e+03
```

## 项目结构

```
├── src/            # 源代码（14 个头文件 + 8 个实现文件）
├── tests/
│   ├── fixtures/   # 测试用例 JSON 文件
│   └── *.cpp       # 7 个测试文件，24 个测试用例
└── docs/superpowers/
    ├── specs/      # 设计规格
    └── plans/      # 实现计划
```

## 许可

MIT
