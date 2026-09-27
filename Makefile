.PHONY: help install test selftest lint clean

CFG ?= configs/minqin.yaml
PY ?= python

help:
	@echo "make install   安装依赖"
	@echo "make test      运行单元测试"
	@echo "make selftest  用合成数据跑通全部 8 个生产脚本（无需真实数据）"
	@echo "make lint      代码检查"
	@echo "make clean     清理缓存"

install:
	$(PY) -m pip install -r requirements.txt

test:
	$(PY) -m pytest

selftest:
	@for f in scripts/0*.py; do \
		echo "=== $$f ==="; \
		$(PY) $$f --config $(CFG) --selftest || exit 1; \
	done

lint:
	$(PY) -m ruff check src scripts tests

clean:
	rm -rf .pytest_cache .ruff_cache
	find . -name __pycache__ -type d -exec rm -rf {} +
