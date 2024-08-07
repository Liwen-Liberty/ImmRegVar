import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from args import Config

config = Config()

class Model(nn.Module):
    def __init__(self, config):
        super(Model, self).__init__()
        self.config = config
        self.embedding_pretrained = config.embedding_pretrained
        
        # Use pre-trained embeddings or initialize embeddings randomly
        if self.embedding_pretrained is not None:
            self.embedding = nn.Embedding.from_pretrained(config.embedding_pretrained, freeze=False)
        else:
            self.embedding = nn.Embedding(config.n_vocab, config.embed)

        # Initialize positional encoding
        self.position_embedding = PositionalEncoding(config.embed, config.pad_size, config.dropout, config.device)
        # Stack multiple encoder layers
        self.encoders = nn.ModuleList([Encoder(config.embed, config.num_head, config.hidden, config.dropout) for _ in range(config.num_encoder)])
        # Initialize custom score mutation layer
        self.my_score_mutation = MyScoreMutation(config.score_mutation_path, config.batch_size, config.n_vocab, config.embed, config.device, config.start)
        # Initialize fully connected layers
        self.fc1 = nn.Linear(config.n_vocab * config.embed, 256)
        self.relu1 = nn.ReLU()
        self.dropout1 = nn.Dropout(p=0.1)
        self.fc2 = nn.Linear(256, config.num_classes)

    def forward(self, x):
        # Generate gene ID matrix
        gene_ids = torch.arange(0, self.config.n_vocab).unsqueeze(0).repeat(x.size(0), 1).to(self.config.device)
        
        # Embedding layer processes the gene ID matrix
        out = self.embedding(gene_ids)
        
        # Positional encoding
        out = self.position_embedding(out)
        
        # Encoder layers processing
        for encoder in self.encoders:
            out, attn_weight = encoder(out)
        
        # Custom score mutation layer processing
        out_emb = self.my_score_mutation(out, x)
        
        # Flatten encoder output
        out = out_emb.view(out_emb.size(0), -1)
        
        # Fully connected layers
        out = self.fc1(out)
        out = self.relu1(out)
        out = self.dropout1(out)
        out = self.fc2(out)
        
        return out, out_emb, attn_weight
    

class MyScoreMutation(nn.Module):
    def __init__(self, score_mutation_path, batch_size, n_vocab, embed, device, start):
        super(MyScoreMutation, self).__init__()
        self.score_mutation_path = score_mutation_path
        self.batch_size = batch_size
        self.n_vocab = n_vocab
        self.embed = embed
        self.device = device
        self.start = start

    def forward(self, out, x):
        final = []
        for i in range(x.shape[0]):
            tem_score_mutation = x[i]
            pair = [1 if score == 1 else 0 for score in tem_score_mutation for _ in range(self.embed)]
            final.append(pair)

        final_tensor = torch.FloatTensor(final).to(self.device)
        out = out.squeeze(-1)
        out = final_tensor.mul(out).to(self.device)

        return out    
    


class Encoder(nn.Module):
    def __init__(self, embed_size, num_head, hidden, dropout):
        super(Encoder, self).__init__()
        # Initialize multi-head attention layer
        self.attention = MultiHeadAttention(embed_size, num_head, dropout)
        # Initialize position-wise feed forward layer
        self.feed_forward = PositionWiseFeedForward(embed_size, hidden, dropout)

    def forward(self, x):
        # Pass through multi-head attention layer
        out, attn_weight= self.attention(x, x, x)
        # Pass through position-wise feed forward layer
        out = self.feed_forward(out)
        return out, attn_weight

class PositionalEncoding(nn.Module):
    def __init__(self, embed, pad_size, dropout, device):
        super(PositionalEncoding, self).__init__()
        self.device = device
        # Generate positional encoding matrix
        pe = torch.zeros(pad_size, embed, device=device)
        position = torch.arange(0, pad_size, dtype=torch.float).unsqueeze(1).to(device)
        div_term = torch.exp(torch.arange(0, embed, 2).float() * (-torch.log(torch.tensor(10000.0)) / embed)).to(device)
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        self.register_buffer('pe', pe)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        # Add positional encoding to input
        x = x + self.pe[:x.size(0), :]
        return self.dropout(x)

class MultiHeadAttention(nn.Module):
    def __init__(self, embed_size, heads, dropout=0.0):
        super(MultiHeadAttention, self).__init__()
        self.embed_size = embed_size
        self.heads = heads
        self.head_dim = embed_size // heads

        assert (
            self.head_dim * heads == embed_size
        ), "Embedding size needs to be divisible by heads"

        # Define linear layers for values, keys, and queries
        self.values = nn.Linear(self.head_dim, embed_size, bias=False)
        self.keys = nn.Linear(self.head_dim, embed_size, bias=False)
        self.queries = nn.Linear(self.head_dim, embed_size, bias=False)
        self.fc_out = nn.Linear(embed_size, embed_size)
        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(embed_size)

    def forward(self, values, keys, query, mask=None):
        N = query.shape[0]
        value_len, key_len, query_len = values.shape[1], keys.shape[1], query.shape[1]

        # Reshape the input to have multiple heads
        values = values.reshape(N, value_len, self.heads, self.head_dim)
        keys = keys.reshape(N, key_len, self.heads, self.head_dim)
        queries = query.reshape(N, query_len, self.heads, self.head_dim)

        # Pass the input through the linear layers
        values = self.values(values)
        keys = self.keys(keys)
        queries = self.queries(queries)

        # Calculate the attention scores
        energy = torch.einsum("nqhd,nkhd->nhqk", [queries, keys])

        if mask is not None:
            energy = energy.masked_fill(mask == 0, float("-1e20"))

        # Apply the softmax function to get the attention weights
        attention = torch.softmax(energy / (self.embed_size ** 0.5), dim=3)

        # Multiply the attention weights with the values
        out = torch.einsum("nhql,nlhd->nqhd", [attention, values]).reshape(
            N, query_len, self.embed_size
        )

        # Pass the result through a linear layer and apply dropout
        out = self.fc_out(out)
        out = self.dropout(out)
        # Apply residual connection and layer normalization
        out = out + query
        out = self.layer_norm(out)
        attention_squeezed = attention.squeeze(1)
        return out, attention_squeezed

class PositionWiseFeedForward(nn.Module):
    def __init__(self, embed_size, hidden, dropout=0.0):
        super(PositionWiseFeedForward, self).__init__()
        # First linear transformation layer of the feed-forward network
        self.fc1 = nn.Linear(embed_size, hidden)
        # Second linear transformation layer of the feed-forward network
        self.fc2 = nn.Linear(hidden, embed_size)
        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(embed_size)

    def forward(self, x):
        # Pass the input through the first linear layer and apply ReLU activation
        out = F.relu(self.fc1(x))
        # Pass the result through the second linear layer and apply dropout
        out = self.fc2(out)
        out = self.dropout(out)
        # Apply residual connection and layer normalization
        out = out + x
        out = self.layer_norm(out)
        return out
