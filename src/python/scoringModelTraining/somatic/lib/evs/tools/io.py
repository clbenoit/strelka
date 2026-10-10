#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

import os
import sys
import pickle
import json


def dict_for_tree(tree):
    """ Make dictionary for tree """
    tree_dict = {}
    for i in range(tree.node_count):
        tree_dict[i] = [int(tree.children_left[i]), int(tree.children_right[i])]
    node_features = tree.feature
    node_thresholds = tree.threshold
    node_count = tree.node_count

    decision_dict = {}
    for i in range(node_count):
        if tree_dict[i][0] == -1 :
            node_features[i] = -1
            node_thresholds[i] = -1

        decision_dict[i] = [int(node_features[i]), float(node_thresholds[i])]
    node_values = tree.value

    # sklearn < 0.22 lays out a classifier's tree.value as (n_nodes, 1, n_classes);
    # sklearn >= 0.22 uses (n_nodes, n_classes). Normalize to the per-node
    # (n_nodes, n_classes) layout so that node_votes[node] = [count(class0),
    # count(class1)] exactly as consumed by the C++ RandomForestModel
    # (vote.left / (vote.left + vote.right)).
    if node_values.ndim == 3:
        node_values = node_values[:, 0, :]
    node_values = node_values.tolist()
    node_values_dict = {}
    for i in range(len(node_values)):
        node_values_dict[i] = node_values[i]
    all_data = {}

    all_data['tree'] = tree_dict
    all_data['decisions'] = decision_dict
    all_data['node_votes'] = node_values_dict

    return all_data


def read_pickled_classifier(pickle_fn):
    """ Read classifier from pickle file """
    with open(pickle_fn, 'rb') as fid:
        clf = pickle.load(fid)
    return clf


def write_classifier_pickle(clf, pickle_fn):
    """ store classifier in pickle file """
    with open(pickle_fn, 'wb') as fid:
        pickle.dump(clf, fid)


def classifier_to_dict(clf):
    """ Convert random forest classifier to  """
    all_trees = []
    for i in range(len(clf.estimators_)):
        all_trees.append(dict_for_tree(clf.estimators_[i].tree_))
    return all_trees


def write_classifier_json(clf, fn):
    """ Write classifier as JSON file """
    all_trees = classifier_to_dict(clf)
    with open(fn, 'w') as fid:
        json.dump(all_trees, fid)
