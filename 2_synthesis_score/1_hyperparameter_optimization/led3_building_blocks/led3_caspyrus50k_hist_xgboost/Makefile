WANDB_ID ?= $(shell bash -c 'read -p "Input Wandb Id: " wandb_id; echo $$wandb_id')
NUMBER_OF_RUNS = 50
ENTITY = <user>
PROJECT = led3_caspyrus50k_hist_xgboost

init_sweep:
	wandb sweep --project $(PROJECT) wandb_config.yaml

run:
	wandb agent --count $(NUMBER_OF_RUNS) $(ENTITY)/$(PROJECT)/$(WANDB_ID)