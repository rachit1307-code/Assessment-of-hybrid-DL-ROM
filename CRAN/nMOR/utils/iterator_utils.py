"""For loading data into NMT models."""
from __future__ import print_function

import collections

import tensorflow as tf

__all__ = ["BatchedInput", "get_iterator", "get_infer_iterator"]


class BatchedInput(
    collections.namedtuple("BatchedInput",
                           ("initializer", "source", "target_output",
                            "sequence_length"))):
  pass


def get_iterator(dataset,
                 batch_size,
                 random_seed,
                 src_max_len=None,
                 tgt_max_len=None,
                 num_parallel_calls=4,
                 output_buffer_size=None):
  if not output_buffer_size:
    output_buffer_size = batch_size * 1000

  # Create source, target, and sequence length dataset
  # TODO: pass input sequence length
  src_tgt_dataset = dataset.map(
      lambda src: (src[0], src, tf.shape(src)[0]),
      num_parallel_calls=num_parallel_calls).prefetch(output_buffer_size)

  # Shuffle around sequence examples
  src_tgt_dataset = src_tgt_dataset.shuffle(output_buffer_size, random_seed)

  # Batch the dataset
  batched_dataset = src_tgt_dataset.batch(batch_size)

  # Make iterator
  batched_iter = batched_dataset.make_initializable_iterator()

  # Get next src, tgt, and seq_len
  src, tgt, seq_len = batched_iter.get_next()

  return BatchedInput(
      initializer=batched_iter.initializer,
      source=src,
      target_output=tgt,
      sequence_length=seq_len)


def get_infer_iterator(dataset, batch_size, num_infer_steps):
  # TODO: pass input sequence length
  dataset = dataset.map(
      lambda src: (src[0], tf.cast(num_infer_steps, tf.int32)))

  batched_dataset = dataset.batch(batch_size)

  batched_iter = batched_dataset.make_initializable_iterator()

  src, seq_len = batched_iter.get_next()

  return BatchedInput(
      initializer=batched_iter.initializer,
      source=src,
      target_output=None,
      sequence_length=seq_len)
