#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
# import torch modules
import torch
import torch.nn as nn
from torch.autograd import grad as torch_grad
from torch.autograd import Variable
# from torch_geometric.nn import GCNConv
# import torch_geometric.nn as gnn
# from torch_geometric.nn import global_mean_pool
# import loss functions
from .losses import *
import numpy as np
import pandas as pd
from scipy.special import softmax
import os
import json
from datetime import date


# refer to: https://github.com/txWang/MOGONET/blob/main/models.py
def xavier_init(m, seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    if type(m) == nn.Linear:
        # print("initial weights")
        # nn.init.xavier_normal_(m.weight)
        nn.init.xavier_uniform_(m.weight)
        if m.bias is not None:
           m.bias.data.fill_(0.0)

########################################
######## Define Layers #################
########################################

# This file refer to https://bytepawn.com/training-a-pytorch-wasserstain-mnist-gan-on-google-colab.html which is code for wasserstein gan, 
# we use it for transfer learning
# https://github.com/EmilienDupont/wgan-gp/blob/master/training.py

def define_linear(input_shape, output_shape, dropout_rate = 0.5, activation = nn.Sigmoid(), seed = 42): 
    if (activation is None):
        layer = nn.Linear(input_shape, output_shape, bias = True)
        
    else:
        if dropout_rate > 0.0:
            layer = nn.Sequential(
            nn.Dropout(dropout_rate),
            nn.Linear(input_shape, output_shape, bias = True),
            nn.BatchNorm1d(output_shape), 
            activation,
            )
        else:
            layer =  nn.Sequential(
            nn.Linear(input_shape, output_shape, bias = True),
            activation)
    layer.apply(lambda x: xavier_init(x, seed))
    return layer
        
# def define_graph_layer(input_shape, output_shape, graph_type  = "GCN", dropout_rate = 0.5, activation = nn.ReLU(inplace = True)):
#     """
#     graph_type is the type of graph layer, GCN as default
#     """
#     if (activation is None):
#         if graph_type == "GCN":
#             layer = GCNConv(input_shape, output_shape)
#         else:
#             raise ValueError("Not Implemented Yet!")
#     else:
#         if graph_type == "GCN":
#             layer = gnn.Sequential('x, edge_index', [
#                 (GCNConv(input_shape, output_shape), 'x, edge_index -> x'),
#                 (gnn.BatchNorm(output_shape), "x -> x"),
#                 (nn.Dropout(dropout_rate), "x -> x"),
#                 (nn.Sigmoid(), "x -> x"),
#         ])
#         else:
#             raise ValueError("Not Implemented Yet!")
#     layer.apply(xavier_init)
#     return layer

class FeatureExtractor(nn.Module): 
    def __init__(self, opt):
        super(FeatureExtractor, self).__init__()
        
        # specify the network structure
        self.emb1_layer = define_linear(opt["input_shape"], opt["feature_dim"], opt["dropout_rate"], nn.Sigmoid(), opt["seed"])
        self.emb2_layer = define_linear(opt["feature_dim"] + opt["input_shape"], opt["feature_dim"], opt["dropout_rate"], nn.Sigmoid(), opt["seed"] + 1)
        self.emb3_layer = define_linear(opt["feature_dim"] * 2 + opt["input_shape"], opt["feature_dim"], opt["dropout_rate"], nn.Sigmoid(), opt["seed"] + 2)
        
        
    def forward(self, input):
        emb1 = self.emb1_layer(input)
        emb2 = self.emb2_layer(torch.cat((input, emb1), axis = 1))
        emb3 = self.emb3_layer(torch.cat((input, emb1, emb2), axis = 1))
        return emb2, emb3


# class GraphFeatureExtractor(nn.Module):
#     def __init__(self, opt):
#         super(GraphFeatureExtractor, self).__init__()
#         # specify the network structure
#         self.emb1_layer = define_graph_layer(opt["input_shape"], opt["feature_dim"], opt["graph_type"], opt["dropout_rate"], nn.Sigmoid())
#         self.emb2_layer = define_graph_layer(opt["feature_dim"] + opt["input_shape"], opt["feature_dim"],  opt["graph_type"], opt["dropout_rate"], nn.Sigmoid())
#         self.emb3_layer = define_graph_layer(opt["feature_dim"] * 2 + opt["input_shape"], opt["feature_dim"], opt["graph_type"], opt["dropout_rate"], nn.Sigmoid())
        

#     def forward(self, x, edge_index, batch):
#         emb1 = self.emb1_layer(x, edge_index)
#         emb2 = self.emb2_layer(torch.cat((x, emb1), axis = 1), edge_index)
#         emb3 = self.emb3_layer(torch.cat((x, emb1, emb2), axis = 1), edge_index)
        
#         emb2 = global_mean_pool(emb2, batch)
#         emb3 = global_mean_pool(emb3, batch)
#         return emb2, emb3
    
class Discriminator(nn.Module):
    def __init__(self, opt):
        super(Discriminator, self).__init__()
        self.model = nn.Sequential(
            *[define_linear(opt["feature_dim"], opt["feature_dim"] // 2, opt["dropout_rate"], nn.LeakyReLU(0.2, inplace = False), opt["seed"]), # will make the code run a little bit slower but don't modify the input
            define_linear(opt["feature_dim"] // 2, opt["feature_dim"] // 2, opt["dropout_rate"], nn.LeakyReLU(0.2, inplace = False), opt["seed"] + 1),
            define_linear(opt["feature_dim"] // 2, 1, 0.0, None, opt["seed"] + 2)]
            )
    def forward(self, input):
        return self.model(input)
    
def create_save_dir(opt):
    if "save_dir" in opt.keys():
        save_dir = opt["save_dir"]
    else:
        save_dir = "checkpoints/{}_{}_{}_{}_{}".format(opt["data_name"], opt["transfer_type"], opt["loss_type"], opt["graph_type"], date.today())
    save_dir = os.path.join(save_dir, "fold_{}_random_seed_{}".format(opt["fold"], opt["seed"])) # folders for different submodels
    return save_dir

############################################################
################ Define Models #############################
############################################################
class BaseModel():
    def __init__(self, opt):
        # set the object's attributes from the optional dictionary
        for key, value in opt.items():
            setattr(self, key, value)
        self.save_dir = create_save_dir(opt)
        os.makedirs(self.save_dir, exist_ok = True)
        # save all configures in a json file
        file_path = os.path.join(self.save_dir, "configs.json")
        with open(file_path, 'w') as json_file:
            json.dump(opt, json_file)
        print("save the configurations into {}".format(file_path))
    
        # add device
        self.device = torch.device('cuda:{}'.format(opt["gpu_id"])) if torch.cuda.is_available() else torch.device('cpu')  # get device name: CPU or GPU
        print("load models on {}".format(self.device))
        # add criterions
        self.criterion_high_norm = nn.MSELoss().to(self.device)
        self.criterion_low_norm = nn.MSELoss().to(self.device)
  
        
    def save_networks(self, epoch):
        # refer to: https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix/blob/master/models/base_model.py
        for name in self.net_name_list:
            net = getattr(self, name)
            if torch.cuda.is_available():
                torch.save(net.cpu().state_dict(), os.path.join(self.save_dir, "{}_net_{}.pth".format(epoch, name)))
                net.cuda(self.gpu_id)
            else:
                torch.save(net.cpu().state_dict(), os.path.join(self.save_dir, "{}_net_{}.pth".format(epoch, name)))
                
    def load_networks(self, epoch):
        # refer to: https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix/blob/master/models/base_model.py
        for name in self.net_name_list:
            load_path = os.path.join(self.save_dir, "{}_net_{}.pth".format(epoch, name))
            net = getattr(self, name)
            if isinstance(net, torch.nn.DataParallel):
                net = net.module
            state_dict = torch.load(load_path, map_location = str(self.device))
            net.load_state_dict(state_dict)
            
    def set_required_grad(self, net_list, requires_grad = True):
        # this part is refer to https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix/blob/da39a525eb793614807db4330cfe9b2157bbe33a/models/base_model.py#L219
        if not isinstance(net_list, list):
            net_list = [net_list]
        for net in net_list:
            if net is not None:
                for param in net.parameters():
                    param.requires_grad = requires_grad
                    
    def set_evaluate_mode(self):
        for net_name in self.net_name_list:
            net = getattr(self, net_name)
            if net is not None:
                net.eval()
                
    def set_train_mode(self):
        for net_name in self.net_name_list:
            net = getattr(self, net_name)
            if net is not None:
                net.train()


    # add functions in basemodel (Oct, 10)
    # refer to: https://github.com/Zeleni9/pytorch-wgan/blob/master/models/wgan_gradient_penalty.py#L296
    def __gradient_penalty(self, real_feat, fake_feat):
        # notice that the real_feat and fake_feat could have different sample size
        # we take minimum of them
        batch_size = min(real_feat.size()[0], fake_feat.size()[0])
        
        # calculate the interpoloation
        alpha = torch.rand(batch_size, 1).to(self.device) # uniformly dist
        alpha = alpha.expand_as(real_feat[:batch_size, ...]).to(self.device)
        interpolated = (alpha * real_feat[:batch_size, ...] + (1 - alpha) * fake_feat[:batch_size, ...]).to(self.device)
        interpolated = Variable(interpolated, requires_grad = True)
        # calculate probability of interpolated examples
        prob_interpolated = self.discrim_layer(interpolated)
        
        # calculate the gradient
        gradients = torch_grad(outputs = prob_interpolated, inputs = interpolated,
                              grad_outputs = torch.ones(prob_interpolated.size()).to(self.device),
                              create_graph = True, retain_graph = True)[0]
        gradients = gradients.view(batch_size, -1)
        
        # calcuate the gradient penalty
        gradients_norm = gradients.norm(2, dim = 1)
        return ((gradients_norm - 1) ** 2).mean()


    # refer to: https://github.com/caogang/wgan-gp/blob/master/gan_mnist.py
    def backward_side(self):
        # optimize discriminator
        self.set_required_grad(self.discrim_layer, True)
        self.optimizer2.zero_grad()
        
        # should be  -torch.mean(real) + torch.mean(fake)
        real_data = Variable(self.low_reso_emb3.data) # since we want to apply classifiers trained on low reso data, we set it as real data
        fake_data = Variable(self.high_reso_emb3.data)
        self.transfer_loss = - torch.mean(self.discrim_layer(real_data)) + torch.mean(self.discrim_layer(fake_data))
        self.transfer_loss += self.gp_weight * self.__gradient_penalty(self.high_reso_emb3.data, self.low_reso_emb3.data)
        self.transfer_loss.backward()
        self.optimizer2.step()

        self.set_required_grad(self.discrim_layer, False)

    # refer to: https://github.com/caogang/wgan-gp/blob/master/gan_mnist.py
    def backward_main(self):
        self.optimizer1.zero_grad()

        self.loss = self.calculate_losses() # this is a function we need to specify for each model
        if self.transfer_type == "Wasserstein":
            self.transfer_loss = - torch.mean(self.discrim_layer(self.high_reso_emb3))
        elif self.transfer_type == "MMD":
            self.transfer_loss = MMD_loss(self.high_reso_emb3, self.low_reso_emb3)

        self.loss += self.lambda3 * self.transfer_loss
        loss_rec = [float(getattr(self, name)) for name in self.loss_names]
        self.loss_rec = pd.concat([self.loss_rec, pd.DataFrame([dict(zip(self.loss_names, loss_rec))])], axis = 0)

        # backward proprogation
        self.loss.backward()
        self.optimizer1.step()
    
    def optimize_parameters(self, epoch):
        self.forward()
        self.backward_main()
        if self.transfer_type == "Wasserstein" and (epoch % self.n_critic == 0): # here we update discriminator much slower since we don't want it to be too smart
                self.forward()
                self.backward_side()

    def linear_eval(self, data_loader, extract_embs = False):
        assert self.graph_type is None # only implement when GNN is not used
        self.set_evaluate_mode() # fix all parameters

        results_df, embs = [], []

        for i, input in enumerate(data_loader):
            self.input = input["data"].float().to(self.device)
            _, self.emb3 = self.feature_extractor_layer(self.input)
            self.pred = self.low_reso_pred_layer(self.emb3).cpu().detach().numpy().squeeze()
            if len(self.pred.shape) > 1:
                self.pred = softmax(self.pred, axis = 1)[:, 1] # for classification, we only extract score for disease
            if "index" in input.keys(): # for sc, st data
                results_df.append(pd.DataFrame({"index": np.array(input["index"]), "hazard": self.pred}))
            elif "pid" in input.keys(): # for pat data
                results_df.append(pd.DataFrame({"pid": np.array(input["pid"]), "hazard": self.pred}))
            else:
                raise ValueError("Can not find index or pid!")

            if extract_embs:
                embs.append(self.emb3.cpu().detach().numpy().squeeze())

        if len(embs) > 0:
            embs = np.vstack(embs)
        return pd.concat(results_df, ignore_index = True), embs

