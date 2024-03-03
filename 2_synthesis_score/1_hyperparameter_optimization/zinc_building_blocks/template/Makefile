WANDB_ID ?= $(shell bash -c 'read -p "Input Wandb Id: " wandb_id; echo $$wandb_id')
NUMBER_OF_RUNS = 50
ENTITY = <user>
PROJECT = <PROJECT_NAME>

init_sweep:
	wandb sweep --project $(PROJECT) wandb_config.yaml

run:
	wandb agent --count $(NUMBER_OF_RUNS) $(ENTITY)/$(PROJECT)/$(WANDB_ID)