import numpy as np
import torch
from load_data import load_data2
from Transformer_model import Model
from args import Config
from sklearn.metrics import r2_score
import torch.nn as nn
import os
from tqdm import tqdm

# Load configuration and data
config = Config()
train, targets = load_data2(config)
script_dir = os.path.dirname(os.path.abspath(__file__))

# Prepare data tensors for training
x_train = torch.from_numpy(train[:]).to(torch.long).to(config.device)
y_train = torch.from_numpy(targets[:]).to(torch.float).to(config.device)

# Initialize model, optimizer, and loss function
model = Model(config)
model.to(config.device)
optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)
loss_func = nn.MSELoss()
model.train()
print('Start iteration....')

empty = torch.empty(0, config.n_vocab * config.embed).to(config.device)
train_acc_temp = 0.5

# Training loop
for step in tqdm(range(config.num_epochs)):
    attention_weights = []
    print('step=', step + 1)
    for batch in range(0, x_train.shape[0], config.batch_size):
        # Forward pass
        train_pre, out_emb, att_weights = model(x_train[batch:batch + config.batch_size])
        
        # Compute loss
        train_loss = loss_func(train_pre.view(-1), y_train[batch:batch + config.batch_size])

        # Compute R^2 score for evaluation
        train_acc = r2_score(y_train[batch:batch + config.batch_size].data.cpu(), train_pre.view(-1).data.cpu())
        train_loss_print = float(train_loss.data.cpu())

        # Backward and optimize
        optimizer.zero_grad()
        train_loss.backward()
        optimizer.step()

        # Collect intermediate embeddings and attention weights
        empty = torch.cat((empty, out_emb), 0)
        attention_weights.extend(att_weights.detach().cpu().numpy())

    # Print training metrics
    print("train_loss:", train_loss_print)
    print("train_acc:", train_acc)

    # Save model and results if accuracy improves
    if train_acc > train_acc_temp:
        np.savetxt(script_dir + '/Results/label.csv', y_train.detach().cpu().numpy(), fmt='%.10f', delimiter=',')
        attention_weights = np.array(attention_weights)
        np.save(script_dir + '/Results/attention_weights.npy', attention_weights)
        content = f'train_loss: = {train_loss_print}\ntrain_acc: = {train_acc}'
        with open(script_dir + '/Results/results.txt', 'w') as f:
            f.write(content)

        torch.save(model.state_dict(), script_dir + '/Results/model_state_dict.pt')

        finaltensor = torch.mm(model.fc2.weight.data, model.fc1.weight.data)
        np.savetxt(script_dir + '/Results/fc1weight.csv', finaltensor.detach().cpu().numpy(), fmt='%.10f', delimiter=',')

        print('Output saved')
        train_acc_temp = train_acc

    empty = torch.empty(0, config.n_vocab * config.embed).to(config.device)


